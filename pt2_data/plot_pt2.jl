using Plots, DelimitedFiles, Polynomials, Statistics, StatsBase, ColorSchemes, Plots.PlotMeasures, LaTeXStrings, FileIO, LsqFit

function sortdata(N, lowbound, upbound)
	#read in all data
	data_init = readdlm("pt2_11/pt2_25_1.txt", ',', Float64, skipstart=1)
	data_init_2 = readdlm("pt2_22/pt2_25_1.txt", ',', Float64, skipstart=1)

	solid = data_init[:,1:2]
	solid_2 = data_init_2[:,1:2]

	#iterate from 2 to N, the last piece of data that will be read-in
	for i = 2:N
		filenam = string("pt2_11/pt2_25_", i, ".txt")
		filenam_2 = string("pt2_22/pt2_25_", i, ".txt")

		data = readdlm(filenam, ',', Float64, skipstart=1)
		data_2 = readdlm(filenam_2, ',', Float64, skipstart=1)

		solid = hcat(solid, data[:,2])
		solid_2 = hcat(solid_2, data_2[:,2])
	end

	#this data is only the solid fractions, not the time
	dset_sol = solid[:,2:N+1]
	dset_sol_2 = solid_2[:,2:N+1]

	#The quartiles are used to calculate the ribbon around the average line seen in the plots
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

	#Plot font
	pfont = "Computer Modern"
	default(fontfamily=pfont)

	#Plot part A
	solidplt = Plots.plot(solid[:,1], avg_sol, grid=false, ribbon = (diffqntl_low, diffqntl_high), fillalpha = 0.5, label  = "", seriescolor = c1) #ribbon and line plot
	scatter!(solid[:,1][1:15:end], avg_sol[1:15:end], grid=false, label = L"r_0/r^* = 1.1", markershape = :rect, markerstrokewidth = 0.5, markersize = 8, legend=:bottomright, seriescolor = c1) #scatter plot for colorblind
	plot!(solid_2[:,1], avg_sol_2, grid=false, ribbon = (diffqntl_low_2, diffqntl_high_2), fillalpha = 0.5, label = "", seriescolor = c2) #ribbon and line plot 2
	scatter!(solid_2[:,1][1:15:end], avg_sol_2[1:15:end], grid=false, label = L"r_0/r^* = 2.2", markershape = :circle, markerstrokewidth = 0.5, markersize = 8, legend=:bottomright, seriescolor = c2) #scatter plot for colorblind2
	plot!(guidefontsize=gfs, legendfontsize=lfs, tickfontsize=tfs, titlefontsize=titlefs, title = "A", titleloc = :left, left_margin = 5mm)
	xlabel!(solidplt, "Time")
	ylabel!(solidplt, "Volume Fraction")

	#Volume Fraction Distrubutions
	lowDistval = floor(Int, 0.1*200)
	highDistval = floor(Int, 0.5*200)

	# Plot parts B and C
	lowDistplot1 = histogram(dset_sol[lowDistval, :], bins=15, alpha=0.5,
	label = L"r_0/r^* = 1.1", titlefontsize = lfs, grid = false,
	xlabel = "                                 Volume Fraction",
	ylabel = "Count", right_margin = 0.5mm, xticks = [0.008, 0.009, 0.011, 0.012],
	top_margin = 0.0mm, title ="B", title_location = :left, ylim = (0, 160))
	lowDistplot2 = histogram(dset_sol_2[lowDistval, :], bins=15, alpha=0.5,
	label = L"r_0/r^* = 2.2", titlefontsize = lfs, grid = false, left_margin = 0.5mm,
	yaxis=false, yticks = false, ylim = (0, 160), color = 2, xticks = 0.08:0.005:0.09,
	top_margin = 0.0mm)
	lowDistplot = plot(lowDistplot1, lowDistplot2)
	plot!(guidefontsize=gfs, legendfontsize=lfs, tickfontsize=tfs, titlefontsize=titlefs)

	highDistplot1 = histogram(dset_sol[highDistval, :], bins=15, alpha=0.5, titlefontsize = lfs, grid = false, xlabel = "                                     Volume Fraction", ylabel = "Count", right_margin = 0.5mm, left_margin = 9mm, ylim = (0, 80), top_margin = 0.0mm, title="C", title_location = :left, legend = false)
	highDistplot2 = histogram(dset_sol_2[highDistval, :], bins=15, alpha=0.5, titlefontsize = lfs, grid = false, left_margin = 0.5mm, yaxis=false, yticks = false, xticks = 0.8:0.1:1.0, color = 2, ylim = (0, 80), top_margin = 0.0mm, legend = false)
	highDistplot = plot(highDistplot1, highDistplot2)
	plot!(guidefontsize=gfs, legendfontsize=lfs, tickfontsize=tfs, titlefontsize=titlefs)

	#Construct the avrami data for plotting
	adj_raw = solid[:,1]
	adj_raw2 = solid_2[:,1]
	index_adj = findall(x -> x > 0.00, adj_raw)
	index_adj2 = findall(x -> x > 0.00, adj_raw2)
	avramiX = log10.(adj_raw[index_adj])
	avramiX_2 = log10.(adj_raw[index_adj2])
	avramidata = log10.(-1 .* log10.(1 .- dset_sol[index_adj, :]))
	avramidata_2 = log10.(-1 .* log10.(1 .- dset_sol_2[index_adj2, :]))
	avramiavg = mean(avramidata, dims=2)
	avramiavg_2 = mean(avramidata_2, dims=2)
	avramistd = std(avramidata, dims=2)
	avramistd_2 = std(avramidata_2, dims=2)

	sfile1 = [avramiX avramiavg]
	sfile2 = [avramiX_2 avramiavg_2]

	sfile_1 = [solid[:,1] avg_sol]
	sfile_2 = [solid_2[:,1] avg_sol_2]

	#This is the fitting using the y bounds
	indexes = findall(x -> x > lowbound && x < upbound, avramiavg)
	indexes_2 = findall(x -> x > lowbound && x < upbound, avramiavg_2)

	avramiY_filt = avramiavg[indexes]
	avramiY_filt_2 = avramiavg_2[indexes_2]
	avramiX_filt = avramiX[indexes]
	avramiX_filt_2 = avramiX_2[indexes_2]

	indsct = unique(round.(Int64, 10 .^ range(0.0, stop=2.3, length = 15)))
	popfirst!(indsct)

	#Plot part D
	avramiplt = plot(avramiX, avramiavg, grid=false, ribbon= 2*avramistd, fillalpha = 0.5, legend = false, c=c1)
	scatter!(avramiX[indsct], avramiavg[indsct], grid=false, markershape = :rect, markerstrokewidth = 0.5, markersize = 8, legend = false, c = c1) #scatter plot for colorblind
	plot!(avramiX_2, avramiavg_2, grid=false, ribbon= 2*avramistd, fillalpha = 0.5, legend = false, c=c2)
	scatter!(avramiX_2[indsct], avramiavg_2[indsct], grid=false, markershape = :circle, markerstrokewidth = 0.5, markersize = 8, legend = false, c = c2)
	plot!(guidefontsize=gfs, legendfontsize=lfs, tickfontsize=tfs, titlefontsize=titlefs ,title = "D", titleloc = :left)

	range1 = findall(x -> x > lowbound && x < upbound, avramiavg)
	range2 = findall(x -> x > lowbound && x < upbound, avramiavg_2)

	xlabel!(avramiplt, "log(t)")
	ylabel!(avramiplt, "log(-log(1-Y))")

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

	#Print the exact average slope value

	println(string("Histogram of Avrami Slopes r* = 2 μ: ", slopeavg, " σ: ", slopestd))
	println(string("Histogram of Avrami Slopes r* = 1 μ: ", slopeavg_2, " σ: ", slopestd_2))

	#Plot part E
	histplt = histogram(hist_slope, bins=15, alpha=0.5, legend = false)
	histogram!(hist_slope_2, bins=15, alpha = 0.5, grid = false, legend = false, title = "E", titleloc = :left)
	histogram!(guidefontsize=gfs, legendfontsize=lfs, tickfontsize=tfs, titlefontsize=titlefs)
	plot!([2.0, 2.0, 2.0], [0, 40, 80], linestyle = :dash, lw = 6, c = :black, legend = false, ylabel = "Count")
	annotate!(2.0, 85, ("Theory", lfs, :black, :center))

	xlabel!("Avrami Constant")

	#Create an array for the reading of various values of phi domains
	curv = fill(0.0, (1251, 1251))
	for b = 0:40:200
	    curv1 = readdlm(string("phi_domains/pt2_25_1_curv_", b, ".txt"), Float64)
	    curv += curv1
	end

	#Plot part F
	curvmin = -0.2
	curvmax = 0.2
	cmap = cgrad([:blue, ColorSchemes.tempo.colors[1], :red], [0.0, 0.45, 0.55, 1.0])
	curvgraph = heatmap(curv, clim=(curvmin, curvmax), seriescolor = cmap, colorbar_title = "Curvature", axis=([], false), title = "F",
			titleloc = :left, layout = grid(1, 2, widths=[0.6, 0.4]))
	heatmap!(titlefontsize=titlefs, colorbar = false, subplot = 1)

	#In order to make the image perfect, I had to do this. Not ideal at all I agree
	img = FileIO.load("cbar.png")
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

	#plot all parts together into one
	plot(solidplt, lowDistplot, highDistplot, avramiplt, histplt, curvgraph, layout = (3,2), size=(1300, 1500))#, background_color=:transparent)
	savefig("../2dp2.png")
end

# 300 data points, lower bound and upper bound are the values of the fitting range in the y-axis
sortdata(300, -2.0, 0.0)
