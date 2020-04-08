γvec = Vector{Float64}()
for u in u_buses
    push!(γvec, γs)
end


fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.055,0.19,0.94,0.795])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$u\$")
#ylim(bottom = 0.015, top = 20)
#ylabel("\$\\chi^{+}_{u}\$")

plot(string(u_buses), χs/γs, color = "lightgreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\beta_{u}\$")
plot(string(u_buses), σ_vec, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\sigma_{u}\$")
plot(string(u_buses), χs, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\chi_{u}\$")
plot(string(u_buses), γvec, color = "cornflowerblue", lw = 1.2, ls = "dashed", label = "\$\\chi\$")

legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9)
savefig(string("plots_final//symmetric.pdf"), format = :pdf)
