using Plots, DelimitedFiles, Polynomials, Statistics, StatsBase, ColorSchemes, Plots.PlotMeasures, LaTeXStrings, FileIO

function sortdata(N, lowbound, upbound)
	#Read in all the data
	data_init = readdlm("pt3_11/pt3_100_1.txt", ',', Float64, skipstart=1)
	data_init_2 = readdlm("pt3_22/pt3_100_1.txt", ',', Float64, skipstart=1)

	solid = data_init[:,1:2]
	solid_2 = data_init_2[:,1:2]

	for i = 2:N
		filenam = string("pt3_11/pt3_100_", i, ".txt")
		filenam_2 = string("pt3_22/pt3_100_", i, ".txt")

		data = readdlm(filenam, ',', Float64, skipstart=1)
		data_2 = readdlm(filenam_2, ',', Float64, skipstart=1)

		solid = hcat(solid, data[:,2])
		solid_2 = hcat(solid_2, data_2[:,2])
	end
	#only volume fraction data contained
	dset_sol = solid[:,2:N+1]
	dset_sol_2 = solid_2[:,2:N+1]

	#Used quartiles for making the ribbons in the plots
	lowQ = 0.1
	highQ = 0.9
	qntl_low = zeros(size(dset_sol)[1])
	qntl_low_2 = zeros(size(dset_sol_2)[1])
	qntl_high = zeros(size(dset_sol)[1])
	qntl_high_2 = zeros(size(dset_sol_2)[1])
	for j = 1:size(dset_sol)[1]
		qntl_low[j] = quantile(dset_sol[j,:], lowQ)
		qntl_low_2[j] = quantile(dset_sol_2[j,:], lowQ)
		qntl_high[j] = quantile(dset_sol[j,:], highQ)
		qntl_high_2[j] = quantile(dset_sol_2[j,:], highQ)
	end

	avg_sol = mean(dset_sol, dims=2)
	avg_sol_2 = mean(dset_sol_2, dims=2)

	std_sol = std(dset_sol, dims=2)
	std_sol_2 = std(dset_sol_2, dims=2)

	diffqntl_low = vec(abs.(avg_sol .- qntl_low))
	diffqntl_high = vec(abs.(qntl_high .- avg_sol))
	diffqntl_low_2 = vec(abs.(avg_sol_2 .- qntl_low_2))
	diffqntl_high_2 = vec(abs.(qntl_high_2 .- avg_sol_2))

	#Font Sizes
	gfs = lfs = tfs = 18
	titlefs = 22

	#Colors
	c1 = 1
	c2 = 2

	#Plot Font
	pfont = "Computer Modern"
	default(fontfamily=pfont, grid=false)

	#Plot part A
	solidplt = plot(solid[:,1], avg_sol, ribbon = (diffqntl_low, diffqntl_high), fillalpha = 0.5, label = "", c = c1)
	scatter!(solid[:,1][1:45:end], avg_sol[1:45:end],label = L"r_0/r^* = 1.1", markershape = :rect, markerstrokewidth = 0.5, markersize = 8, c = c1)
	plot!(solid_2[:,1], avg_sol_2, grid=false, ribbon = (diffqntl_low_2, diffqntl_high_2), fillalpha = 0.5, label = "", c = c2)
	scatter!(solid_2[:,1][1:45:end], avg_sol_2[1:45:end],label = L"r_0/r^* = 2.2", markershape = :circle, markerstrokewidth = 0.5, markersize = 8, c = c2, legend = :bottomright)
	plot!(guidefontsize=gfs, legendfontsize=lfs, tickfontsize=tfs, titlefontsize=titlefs, title = "A", titleloc = :left, left_margin = 5mm)
	xlabel!(solidplt, "Time")
	ylabel!(solidplt, "Volume Fraction")

	# Distrubution of Volume Fraction

	lowDistval = floor(Int, 100)
	highDistval = floor(Int, 250)

	#Plot part B
	lowDistplot = histogram(dset_sol[lowDistval, :], bins=15, alpha=0.5, label=L"r_0/r^* = 1.1")
	histogram!(dset_sol_2[lowDistval, :], bins=15, alpha=0.5, label=L"r_0/r^* = 2.2", title = "B", titleloc = :left, right_margin = 3mm, ylim = (0, 90))
	plot!(guidefontsize=gfs, legendfontsize=lfs, tickfontsize=tfs, titlefontsize=titlefs)
	xlabel!("Volume Fraction")
	ylabel!("Count")

	#Plot part C
	highDistplot1 = histogram(dset_sol[highDistval, :], bins=9, alpha=0.5, titlefontsize = lfs, grid = false, xlabel = "                            Volume Fraction", ylabel = "Count", left_margin = 9.5mm, right_margin = 5.0mm, top_margin = 0.0mm, title = "C", title_location = :left, legend = false)
	highDistplot2 = histogram(dset_sol_2[highDistval, :], bins=9, alpha=0.5, titlefontsize = lfs, grid = false, left_margin = 0.5mm, yaxis= false, yticks = false, color = 2, top_margin = 0.0mm, legend = false)
	highDistplot = plot(highDistplot1, highDistplot2)
	plot!(guidefontsize=gfs, legendfontsize=lfs, tickfontsize=tfs, titlefontsize=titlefs)

	#Create Avrami Data
	avramidata = log10.(-1 .* log10.(1 .- dset_sol))
	avramidata_2 = log10.(-1 .* log10.(1 .- dset_sol_2))

	avramiavg = mean(avramidata, dims=2)
	avramiavg_2 = mean(avramidata_2, dims=2)
	avramistd = std(avramidata, dims=2)
	avramistd_2 = std(avramidata_2, dims=2)
	avramiX = log10.(solid[:,1])
	avramiX_2 = log10.(solid_2[:,1])

	#Index according to the fitting range
	indexes = findall(x -> x > lowbound && x < upbound, avramiavg)
	indexes_2 = findall(x -> x > lowbound && x < upbound, avramiavg_2)

	avramiY_filt = avramiavg[indexes]
	avramiY_filt_2 = avramiavg[indexes_2]
	avramiX_filt = avramiX[indexes]
	avramiX_filt_2 = avramiX[indexes_2]

	#find quantiles
	qavr_low = zeros(size(avramidata)[1])
	qavr_low_2 = zeros(size(avramidata_2)[1])
	qavr_high = zeros(size(avramidata)[1])
	qavr_high_2 = zeros(size(avramidata_2)[1])
	avr = zeros(size(avramidata)[1])
	avr_2 = zeros(size(avramidata_2)[1])
	for j = 1:size(avramidata)[1]
		qavr_low[j] = quantile(avramidata[j,:], lowQ)
		qavr_low_2[j] = quantile(avramidata_2[j,:], lowQ)
		qavr_high[j] = quantile(avramidata[j,:], highQ)
		qavr_high_2[j] = quantile(avramidata_2[j,:], highQ)
		avr[j] = quantile(avramidata[j,:], 0.5)
		avr_2[j] = quantile(avramidata_2[j,:], 0.5)

	end

	diffqavr_low = vec(abs.(avr .- qavr_low))
	diffqavr_high = vec(abs.(qavr_high .- avr))
	diffqavr_low_2 = vec(abs.(avr_2 .- qavr_low_2))
	diffqavr_high_2 = vec(abs.(qavr_high_2 .- avr_2))

	#Sort through data to eliminate NaNs and Infinite values
	idx_avr_low = findall(x -> x == NaN || x == Inf || x == -Inf, diffqavr_low)
	idx_avr_high = findall(x -> x == NaN || x == Inf || x == -Inf, diffqavr_high)
	idx_avr_low_2 = findall(x -> x == NaN || x == Inf || x == -Inf, diffqavr_low_2)
	idx_avr_high_2 = findall(x -> x == NaN || x == Inf || x == -Inf, diffqavr_high_2)
	diffqavr_low[idx_avr_low] .= 1.0
	diffqavr_high[idx_avr_high] .= 0.0
	diffqavr_low_2[idx_avr_low_2] .= 1.0
	diffqavr_high_2[idx_avr_high_2] .= 0.0

	indsct = unique(round.(Int64, 10 .^ range(0.0, stop=2.77, length = 15)))

	#Plot part D
	avramiplt = plot(avramiX, avr, ribbon= (diffqavr_low, diffqavr_high), fillalpha = 0.5, legend=false, c = c1)
	scatter!(avramiX[indsct], avr[indsct], markershape = :rect, markerstrokewidth = 0.5, markersize = 8, legend = false, c = c1)
	plot!(avramiX_2, avr_2, grid=false, ribbon= (diffqavr_low_2, diffqavr_high_2), fillalpha = 0.5, legend=false, c = c2)
	scatter!(avramiX_2[indsct], avr_2[indsct], markershape = :circle, markerstrokewidth = 0.5, markersize = 8, c = c2, title = "D", titleloc = :left, ylims = (-5.05, 1))
	plot!(guidefontsize=gfs, legendfontsize=lfs, tickfontsize=tfs, titlefontsize=titlefs, xlims = (0.6, 3.0), right_margin = 3mm)

	f = Polynomials.fit(avramiX_filt, avramiY_filt, 1)
	f_2 = Polynomials.fit(avramiX_filt_2, avramiY_filt_2, 1)
	slope = f[1]
	slope_2 = f_2[1]
	int = f[0]
	int_2 = f_2[0]

	println("Testing")
	println(string("Graph of Avrami Slopes r*=2 μ: ", slope))
	println(string("Graph of Avrami Slopes r*=1 μ: ", slope_2))

	xlabel!("log(t)")
	ylabel!("log(-log(1-Y))")

	#Histogram of Avrami Slopes
	hist_slope = zeros(N)
	hist_slope_2 = zeros(N)
	for i = 1:N
		avData = avramidata[:,i]
		avData_2 = avramidata_2[:,i]
		indexes = findall(x -> x > lowbound && x < upbound, avData)
		indexes_2 = findall(x -> x > lowbound && x < upbound, avData_2)
		avramiX_fil = avramiX[indexes]
		avramiX_fil_2 = avramiX_2[indexes_2]
		avramiY_fil = avData[indexes]
		avramiY_fil_2 = avData_2[indexes_2]

		fit_iter = Polynomials.fit(avramiX_fil, avramiY_fil, 1)
		fit_iter_2 = Polynomials.fit(avramiX_fil_2, avramiY_fil_2, 1)
		hist_slope[i] = fit_iter[1]
		hist_slope_2[i] = fit_iter_2[1]
	end

	slopeavg = mean(hist_slope)
	slopeavg_2 = mean(hist_slope_2)
    slopestd = std(hist_slope)
	slopestd_2 = std(hist_slope_2)

	#Print exact averages
	println(string("Histogram of Avrami Slopes r* = 2 μ: ", slopeavg, " σ: ", slopestd))
	println(string("Histogram of Avrami Slopes r* = 1 μ: ", slopeavg_2, " σ: ", slopestd_2))

	#Plot part E
	histplt = histogram(hist_slope, bins=15, alpha=0.5, legend = false)
	histogram!(hist_slope_2, bins=30, alpha = 0.5, legend = false, grid = false, title = "E", titleloc = :left, xlims=(2.4,4.0))
	histogram!(guidefontsize=gfs, legendfontsize=lfs, tickfontsize=tfs, titlefontsize=titlefs)
	plot!([3.0, 3.0, 3.0], [0, 40, 58], linestyle = :dash, lw = 6, c = :black, legend = false, ylabel = "Count")
	annotate!(3.0, 62.5, ("Theory", lfs, :black, :center))
	xlabel!("Avrami Constant")

	#Create array for curvature data
	curv = fill(0.0, (2501, 2501))
	for b = 0:50:600
	    curv1 = readdlm(string("phi_domains/pt3_100_1_curv_", b, ".txt"), Float64)
	    curv += curv1
	end

	#Plot part F
	curvmin = -0.2
	curvmax = 0.2
	cmap = cgrad([:blue, ColorSchemes.tempo.colors[1], :red], [0.0, 0.47, 0.53, 1.0])
	curvgraph = heatmap(curv, clim=(curvmin, curvmax), seriescolor = cmap, colorbar_title = "Curvature", axis=([], false),
			title = "F", titleloc = :left, layout = grid(1, 2, widths=[0.6, 0.4]))
	heatmap!(titlefontsize=titlefs, colorbar = false, subplot = 1)

	#Again, not ideal at all but this was necessary to make the graph exactly how I wanted it
	img = FileIO.load("../pt2_data/cbar.png")
	plot!(img, axis = ([], false), subplot =2, left_margin = -35.0mm)
	ann_x = 37
	annotate!(ann_x+20, -2, ("0.20", lfs, :black, :left), subplot = 2)
	annotate!(ann_x, -16, ("_", lfs, :black, :left), subplot = 2)
	annotate!(ann_x+20, 90.5, ("0.10", lfs, :black, :left), subplot = 2)
	annotate!(ann_x, 77, ("_", lfs, :black, :left), subplot = 2)
	annotate!(ann_x+20, 183, ("0", lfs, :black, :left), subplot = 2)
	annotate!(ann_x, 169, ("_", lfs, :black, :left), subplot = 2)
	annotate!(ann_x+20, 280, ("-0.10", lfs, :black, :left), subplot = 2)
	annotate!(ann_x, 262, ("_", lfs, :black, :left), subplot = 2)
	annotate!(ann_x+20, 375, ("-0.20", lfs, :black, :left), subplot = 2)
	annotate!(ann_x, 355, ("_", lfs, :black, :left), subplot = 2)
	annotate!(140, 245, ("Curvature", lfs, :black, :left), annotationrotation = 90.0, subplot = 2)

	#combine all plots into one
	plot(solidplt, lowDistplot, highDistplot, avramiplt, histplt, curvgraph, layout = (3,2), size=(1300, 1500))
	savefig("../2dp3.png")
end

sortdata(300, -2.0, 0.0)
