fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.095,0.19,0.90,0.795])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$i\$")
xlim(left = 1, right = n_generators)
#ylim(bottom = 0.015, top = 20)
ylabel("\$\\delta_{i}\$")

plot([i for i in 1:n_generators], abs(cc1_ab), color = "navy", lw = 0.6, marker = "D", ms = 2.0, markeredgewidth = 0.6, mec = "navy", mfc = "white", ls = "dotted", label = "\$\\overline{\\delta}_{i}\$")
plot([i for i in 1:n_generators], abs(cc1_ab_c), color = "mediumblue", lw = 0.6, marker = "D", ms = 2.0, markeredgewidth = 0.6, mec = "mediumblue", mfc = "white", ls = "dotted", label = "\$\\overline{\\delta}_{i}\$-c")
plot([i for i in 1:n_generators], abs(cc2_ab), color = "teal", lw = 0.6, marker = "D", ms = 2.0, markeredgewidth = 0.6, mec = "teal", mfc = "white", ls = "dotted", label = "\$\\underline{\\delta}_{i}\$")
plot([i for i in 1:n_generators], abs(cc2_ab_c), color = "lightseagreen", lw = 0.6, marker = "D", ms = 2.0, markeredgewidth = 0.6, mec = "lightseagreen", mfc = "white", ls = "dotted", label = "\$\\underline{\\delta}_{i}\$-c")

legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9)
savefig(string("plots_final//delta.pdf"), format = :pdf)
