#Nucleation Benchmark Code
#CUDA Implimentation - Julia
#Author: Jose Mancias

using CUDA, DelimitedFiles

# The primary kernal for calculating an iteration of a phase field calculation
function iterat(phi::CuDeviceMatrix{Float32,1}, norm, curv_2, Nl::Int64, Ml::Int64, xiter::Float64, Titer::Float64, dfl::Float64)
    idx_i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    idx_j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if idx_i < Nl && idx_j < Ml && idx_i > 1 && idx_j > 1

        #Calculate the necessary physics
		curr_phi = phi[idx_i, idx_j]
        g_prime = 2*phi[idx_i,idx_j]*(phi[idx_i,idx_j]-1)*(2*phi[idx_i,idx_j]-1)
        p_prime = 30*phi[idx_i,idx_j]^2*(phi[idx_i,idx_j]^2-2*phi[idx_i,idx_j]+1)

        laplacian_phi = (phi[idx_i+1, idx_j]+phi[idx_i-1, idx_j]+phi[idx_i, idx_j-1]+phi[idx_i, idx_j+1]-4.0*phi[idx_i,idx_j])/(xiter*xiter)


		#curvature calculation
		gradient_x = (phi[idx_i+1, idx_j]-phi[idx_i-1, idx_j])/(2.0*xiter)
        gradient_y = (phi[idx_i, idx_j+1]-phi[idx_i, idx_j-1])/(2.0*xiter)

		n_x = gradient_x / norm[idx_i, idx_j]
		n_y = gradient_y / norm[idx_i, idx_j]

		A_x = (norm[idx_i+1, idx_j]-norm[idx_i-1, idx_j])/(2.0*xiter)
		A_y = (norm[idx_i, idx_j+1]-norm[idx_i, idx_j-1])/(2.0*xiter)

		B = (n_x*A_x)+(n_y*A_y)

		#Don't allow a division by 0
		if(norm[idx_i, idx_j] <= 0.0001)
			curv_2[idx_i, idx_j] = 0.0
		else
			curv_2[idx_i, idx_j] = -1*((laplacian_phi - B) / norm[idx_i, idx_j])*(exp(-1*(phi[idx_i, idx_j]-0.5)^2 / 0.01))
		end

		#Calculate the new phi value
		@inbounds phi[idx_i,idx_j] = curr_phi+Titer*(laplacian_phi - g_prime + dfl*p_prime)
    end

	# Periodic Boundary Conditions
    if idx_i == 2 && idx_j <= Ml
		phi[Nl, idx_j] = phi[2, idx_j]
    end
    if idx_j == 2 && idx_i <= Nl
		phi[idx_i, Ml] = phi[idx_i, 2]
    end
    if idx_i == Nl-1 && idx_j <= Ml
		phi[1, idx_j] = phi[Nl-1, idx_j]
    end
    if idx_j == Ml-1 && idx_i <= Nl
		phi[idx_i, 1]  = phi[idx_i, Ml-1]
    end

    return
end

#This function calculates the norm, a value used in the calculation of the curvature
function normCalc(phi::CuDeviceMatrix{Float32,1}, norm, Nl, Ml, xiter)
    idx_i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    idx_j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

    if idx_i < Nl && idx_j < Ml && idx_i > 1 && idx_j > 1
		gradient_x = (phi[idx_i+1, idx_j]-phi[idx_i-1, idx_j])/(2.0*xiter)
        gradient_y = (phi[idx_i, idx_j+1]-phi[idx_i, idx_j-1])/(2.0*xiter)
        norm[idx_i, idx_j] = sqrt(gradient_x^2 + gradient_y^2)
    end
    return
end

#initialize the seeds from the list of seeds developed
function initSeedsKern(phi::CuDeviceMatrix{Float32,1}, Nl, Ml, xiter, xpos, ypos, rinit)
    idx_i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    idx_j = (blockIdx().y - 1) * blockDim().y + threadIdx().y

	#convert to "dx coordinates"
	x_p = idx_i*xiter
	y_p = idx_j*xiter

    if idx_i <= Nl && idx_j <= Ml && idx_i >= 1 && idx_j >= 1
		#add phi value of a single seed to the phi domain
		phiold = phi[idx_i, idx_j]
		radius = sqrt((x_p - xpos)^2 + (y_p - ypos)^2)
		phi1 = (1.0 - tanh((radius - rinit)/(sqrt(2))))/2
		phi[idx_i, idx_j] = phiold + phi1

		#ensure no value of phi is greater than the limit, 1.0
		if phi[idx_i, idx_j] > 1.0
			phi[idx_i, idx_j] = 1.0
		end
    end
    return
end

#Represents the method for a single run of a benchmark problem starting at t=0 till end.
function single(N::Int64, M::Int64, dx::Float64, Time::Float64, filename::String, seeds, df::Float64)
	#GPU values
	nthreads = 16
    nblocksx = ceil(Int, N/nthreads)
    nblocksy = ceil(Int, M/nthreads)

	#Time values
	dt = round((dx*dx)/16, digits=3) #make dt smaller than the limit of dx*dx / 4
    Tsave = round.(collect(0.0:1.0:Time), digits=6) #What values will the data be recorded at
	TsaveFull = [80, 85, 90] #These values describe when the entire PHI array will be recorded

	#define fields
    phiCurr_h = fill(0.0f0, (N,M)) #host PHI array
    phiCurr = CUDA.fill(0.0f0, (N,M)) #GPU PHI array

    curv_2_h = similar(phiCurr_h) #host curvature array
    curv_2 = similar(phiCurr) #GPU curvature array
	curv_pen = fill(0.0f0, (N-2, M-2)) #

    norm = similar(phiCurr)

	#loop to march in time
	indseed = 1

	io = open(filename, "w")
    write(io, "time, volume fraction, curvature(filtered), pos, neg \n")
    close(io)
    for t = 0.0:dt:Time
		while(indseed <= length(seeds[:,1]) && t-dt < seeds[indseed, 4] && t >= seeds[indseed, 4])

			xpos = seeds[indseed, 1]
			ypos = seeds[indseed, 2]
			r0 = seeds[indseed, 3]
			@cuda blocks=(nblocksx, nblocksy) threads=(nthreads, nthreads) initSeedsKern(phiCurr, N, M, dx, xpos, ypos, r0)

	    	indseed = indseed + 1
		end

		@cuda blocks=(nblocksx, nblocksy) threads=(nthreads, nthreads) normCalc(phiCurr, norm, N, M, dx)

	    @cuda blocks=(nblocksx, nblocksy) threads=(nthreads, nthreads) iterat(phiCurr, norm, curv_2, N, M, dx, dt, df)


		if(t in Tsave)
            Vf = sum(phiCurr[2:(N-1), 2:(M-1)]) / ((N-2)*(M-2)) #Find volume fraction

			# remove NaN values
			curv_pen[isnan.(curv_pen)] .= 0.0
			pos = filter(x->x>0, curv_pen) #Part of a curvature analysis where the postive values of curvature were calculated
			neg = filter(x->x<0, curv_pen) #Part of a curvature analysis again, this time negative curvature
			Cv = sum(curv_pen) / ((N-2)*(M-2)) #Find total curvature present
			Pv = sum(pos) / ((N-2)*(M-2)) #Postive value of curvature
			Nv = sum(neg) / ((N-2)*(M-2)) #Negative value of curvature

	    	fileStream = open(filename, "a")
	    	write(fileStream, string(t, ",", Vf, ",", Cv, ",", Pv, ",", Nv, "\n"))
	    	close(fileStream)

		end

		if(t in TsaveFull)
			t_int = floor(Int, t)

			# copy back to cpu
			copyto!(phiCurr_h, phiCurr)

			#This large number is used to help in the naming of the file
			sizeInd = floor(Int64, seeds[3]*100)
			phi_save = string("phi_", t_int, "_", sizeInd, ".txt")

			open(phi_save, "w") do io
				writedlm(io, phiCurr_h)
			end

		end

	end
	println(indseed)

    return
end

#Part 1 of the benchmark
function main_one()
	#initial conditions
    dx = 0.4
    Lx = Ly = 100
    rstar = 5
	r0 = 5
    df = sqrt(2)/30
    timeEnd = 100.0

	# each run, the single method is called. Single is called using the following parameters:
	#single(Number of mesh points in X, Y, dx, timeEnd, stringForSavingFile, [X position of seeds | Y positions | r0 value | time], df)
    single(floor(Int, Lx/dx+1), floor(Int, Ly/dx+1), dx, timeEnd, "nucPt1_1.txt", [Lx/2 Ly/2 r0 0.0], df)
    single(floor(Int, Lx/dx+1), floor(Int, Ly/dx+1), dx, timeEnd, "nucPt1_101.txt", [Lx/2 Ly/2 r0*1.01 0.0], df)
    single(floor(Int, Lx/dx+1), floor(Int, Ly/dx+1), dx, timeEnd, "nucPt1_099.txt", [Lx/2 Ly/2 r0*0.99 0.0], df)
    single(floor(Int, Lx/dx+1), floor(Int, Ly/dx+1), dx, timeEnd, "nucPt1_11.txt", [Lx/2 Ly/2 r0*1.1 0.0], df)
    single(floor(Int, Lx/dx+1), floor(Int, Ly/dx+1), dx, timeEnd, "nucPt1_22.txt", [Lx/2 Ly/2 r0*2.2 0.0], df)
end

#Part 2 of the benchmark problem
function main_two()
	#initial conditions
    dx = 0.4
    num_sol = 300 #Number of runs
    Lx = Ly = 500
    rstar = 2
    df = sqrt(2)/6
    timeEnd = 200.0
    r0 = 2.2
    seedsnum = 25

    for i = 1:num_sol
		seedtimes = zeros(seedsnum,1) #an array of 0s for the time portion of the seeds
		seedradii = ones(seedsnum,1) .* r0 #an array of the value of seed radius
		seedposition = rand(Float64, (seedsnum,2)).*Lx #a randomized array for calculating 2 columns of random values between 0 and Lx
        seeds = hcat(seedposition, seedradii, seedtimes) #put the three arrays above together
		savepath = string("pt2_", seedsnum, "_", i, ".txt")
		single(floor(Int, Lx/dx+1), floor(Int, Ly/dx+1), dx, timeEnd, savepath, seeds, df)
    end
end

#Part 3 of the benchmark problem
function main_three()
	#Initial conditions
   dx = 0.4
   num_sol = 300
   Lx = Ly = 1000
   rstar = 1
   df = sqrt(2)/6
   timeEnd = 600.0
   r0 = 2.2
   seedsnum = 100

   for i = 1:num_sol
	   seedtimes = rand(Float64, (seedsnum,1)).*timeEnd #rand times in part 3
	   seedtimes = sort(seedtimes, dims=1) #sort the above array from lowest times to highest
	   seedradii = ones(seedsnum,1) .* r0
	   seedposition = rand(Float64, (seedsnum,2)).*Lx
	   seeds = hcat(seedposition, seedradii, seedtimes)
	   savepath = string("pt3_", seedsnum, "_", i, ".txt")
	   single(floor(Int, Lx/dx+1), floor(Int, Ly/dx+1), dx, timeEnd, savepath, seeds, df)
   end

end

#@time means that each function will be timed and you will know how long each part takes
@time main_one()
@time main_two()
@time main_three()
