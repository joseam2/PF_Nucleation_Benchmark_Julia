using Plots, DelimitedFiles, Plots.PlotMeasures, ColorSchemes, FileIO
#Graph 1 in Paper

v099 = readdlm("pt1_099.txt", ',', Float64, skipstart=1)
v100 = readdlm("pt1_100.txt", ',', Float64, skipstart=1)
v101 = readdlm("pt1_101.txt", ',', Float64, skipstart=1)
v11 = readdlm("pt1_11.txt", ',', Float64, skipstart=1)
v200 = readdlm("nucPt1_2.txt", ',', Float64, skipstart=1)

#Plot font
pfont = "Computer Modern"
ft = 14
default(fontfamily=pfont, guidefontsize = ft, legendfontsize = ft, tickfontsize = ft)

solid = plot(v099[:,1], v099[:,2], label = "", ls = :solid, linewidth = 4, c = 1, markershape = :none)
scatter!(v099[:,1][1:15:end], v099[:,2][1:15:end], label = "0.99r*", markershape = :rect, markerstrokewidth = 0.5, markersize = 8, c = 1)
plot!(v100[:,1], v100[:,2], label = "", ls = :solid, linewidth = 4, c = 2, markershape = :none)
scatter!(v100[:,1][1:15:end], v100[:,2][1:15:end], label = "1.00r*", markershape = :circle, markerstrokewidth = 0.5, markersize = 8, c = 2)
plot!(v101[:,1], v101[:,2], label = "", ls = :solid, linewidth = 4, c = 3, markershape = :none)
scatter!(v101[:,1][1:15:end], v101[:,2][1:15:end], label = "1.01r*", markershape = :dtriangle, markerstrokewidth = 0.5, markersize = 8, c = 3)
plot!(v11[:,1], v11[:,2], label = "", ls = :solid, linewidth = 4, legend = :topright, c = 4, markershape = :none)
scatter!(v11[:,1][1:15:end], v11[:,2][1:15:end], label = "1.1r*", markershape = :utriangle, markerstrokewidth = 0.5, markersize = 8, c = 4)
plot!(v200[:,1], v200[:,2], label = "", ls = :solid, linewidth = 4, c = 6, markershape = :none)
scatter!(v200[:,1][1:15:end], v200[:,2][1:15:end], label = "2.0r*", markershape = :diamond, markerstrokewidth = 0.5,
        markersize = 8, c = 6, grid = false, legend_hfactor=1.25, extra_kwargs=:subplot, size = (850, 600), left_margin = 5mm, ylims = (0, 0.073))

xlabel!("Time")
ylabel!("Volume Fraction")

#Add domain image
img = FileIO.load("p1_domain.png")
plot!(img, axis = ([], false), annotate = (50, 80, ("t = 90", ft+1, :black, :left)),
                inset = (1, bbox(0.15, 0.032, 0.9, 0.45)), subplot =2)
savefig("../part1_volfrac.png")
