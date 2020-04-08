γvec = Vector{Float64}()
γmvec = Vector{Float64}()
γpvec = Vector{Float64}()
for u in u_buses
    push!(γvec, γs)
    push!(γpvec, γp)
    push!(γmvec, γm)
end

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.095,0.19,0.9,0.795])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$u\$")
#ylim(bottom = 0.015, top = 20)
ylabel("\$\\chi^{\\pm}_{(u)},~\\beta^{\\pm}\$")

#plot(string(u_buses), abs(γvec), color = "cornflowerblue", lw = 1.2, ls = "solid", label = "\$\\chi\$")
plot(string(u_buses), abs(γpvec), color = "blue", lw = 0.8, marker = "D", ms = 3.0, markeredgewidth = 0.8, mec = "blue", mfc = "white", ls = "dashed", label = "\$\\chi^{\\pm}\$")
#plot(string(u_buses), abs(γmvec), color = "cornflowerblue", lw = 1.2, ls = "dashed", label = "\$\\chi^{-}\$")
plot(string(u_buses), abs(χp), color = "cornflowerblue", lw = 0.8, ls = "dashed", marker = "D", markeredgewidth = 0.8, ms = 3.0, mfc = "white", label = "\$\\chi^{\\pm}_{u}\$")
plot(string(u_buses), abs(χm/γm), color = "lightseagreen", lw = 0.8, ls = "dashed", marker = "D", markeredgewidth = 0.8, ms = 3.0, mfc = "white", label = "\$\\beta^{\\pm}_{u}\$")
#plot(string(u_buses), χm/γm, color = "lightgreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\beta^{-}_{u}\$")
#plot(string(u_buses), abs(χs/γs), color = "lightgreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\beta_{u}\$")
plot(string(u_buses), σm, color = "lightgreen", lw = 0.8, ls = "dashed", marker = "D", markeredgewidth = 0.8, ms = 3.0, mfc = "white", label = "\$\\sigma_{u}\$")

#plot(string(u_buses), abs(χm), color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\chi^{-}_{u}\$")
#plot(string(u_buses), abs(χs), color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\chi_{u}\$")

legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9)
savefig(string("plots_final//symmetric_asymmetric.pdf"), format = :pdf)
