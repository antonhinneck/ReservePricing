using Colors
using PyPlot

fig = figure(figsize=(1.6, 1.5))
rc("font", family = "serif", style = "italic", size = 11)
rc("text", usetex = true)
rc("lines", linewidth = 1)
rc("text.latex", preamble = "\\usepackage{xcolor}")

ax = fig.add_axes([0.18,0.14,0.8,0.825])
# grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ylabel("\$[\\frac{\\\$}{h}]\$")
yticks([])

x = [1, 2, 3, 4]
colors = ["lightseagreen", "seagreen", "navy", "blue"]
labls = ["\$\\hat{\\chi}\$", "\$\\chi\$", "\$\\chi^{+}\$", "\$\\chi^{-}\$"]
xticks(x, labls)
ylim(bottom = 0.0, top = 80.0)

#ax.get_yaxis().set_visible(false)

bar(x, [χ0, χ1, χm, χp], align = "center", color = colors, width = 0.5)

for (i, val) in enumerate([χ0, χ1, χm, χp])
    annotate(string(round(val, digits = 1)), (i - 0.42, val + 2.0))
end

savefig(string("plots_final//chi_sw.pdf"), format = :pdf)

#-----------------
#-----------------

fig = figure(figsize=(5.0, 1.5))
rc("font", family = "serif", style = "italic", size = 11)
rc("text", usetex = true)
rc("lines", linewidth = 1)
rc("text.latex", preamble = "\\usepackage{xcolor}")

ax = fig.add_axes([0.1,0.24,0.89,0.74])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ylabel("\$[\\frac{\\\$}{h}]\$")
xlabel("\$u\$")

#ylim(bottom = -1.0, top = 25)

#ax.set_yscale("log")
colors = ["lightseagreen", "seagreen", "navy", "blue"]
# labls = ["\$\\hat{\\chi}\$", "\$\\chi\$", "\$\\chi^{+}\$", "\$\\chi^{-}\$"]
# ylim(bottom = 0.0, top = 80.0)
# ax.get_yaxis().set_visible(false)

x = [i for i in 1:length(uRESs)]
labls = [string(i) for i in 1:length(uRESs)]
xticks(x, labls)

plot(x, χ0u, color = colors[1], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\bar{\\chi}_{u}\$")
plot(x, χ1u, color = colors[2], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\chi_{u}\$")
# plot(x, μ, color = colors[3], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\mu\$")
# plot(x, σ, color = colors[4], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\mu\$")
# plot(x, χ1up, color = colors[3], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\chi^{+}_{u}\$")
# plot(x, χ1um, color = colors[4], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\chi^{-}_{u}\$")

# for (i, val) in enumerate([χ0, χ1, χm, χp])
#     annotate(string(round(val, digits = 1)), (i - 0.4, val + 2.0))
# end

legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 4)
savefig(string("plots_final//chi_n2n.pdf"), format = :pdf)

#-------------
#-------------

fig = figure(figsize=(5.0, 1.5))
rc("font", family = "serif", style = "italic", size = 11)
rc("text", usetex = true)
rc("lines", linewidth = 1)
rc("text.latex", preamble = "\\usepackage{xcolor}")

ax = fig.add_axes([0.105,0.24,0.89,0.74])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ylabel("\$[\\frac{\\\$}{h}]\$")
xlabel("\$u\$")

#ylim(bottom = -1.0, top = 25)

#ax.set_yscale("log")
colors = ["lightseagreen", "seagreen", "navy", "blue"]
# labls = ["\$\\hat{\\chi}\$", "\$\\chi\$", "\$\\chi^{+}\$", "\$\\chi^{-}\$"]
# ylim(bottom = 0.0, top = 80.0)
# ax.get_yaxis().set_visible(false)

x = [i for i in 1:length(uRESs)]
labls = [string(i) for i in 1:length(uRESs)]
xticks(x, labls)

plot(x, μ, color = colors[1], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\mu\$")
plot(x, σ, color = colors[2], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\sigma\$")
# plot(x, χ1up, color = colors[3], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\chi^{+}_{u}\$")
# plot(x, χ1um, color = colors[4], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\chi^{-}_{u}\$")

# for (i, val) in enumerate([χ0, χ1, χm, χp])
#     annotate(string(round(val, digits = 1)), (i - 0.4, val + 2.0))
# end

legend(loc = "upper left", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 4)
savefig(string("plots_final//params.pdf"), format = :pdf)

#-----------------
#-----------------

fig = figure(figsize=(5.0, 1.5))
rc("font", family = "serif", style = "italic", size = 11)
rc("text", usetex = true)
rc("lines", linewidth = 1)
rc("text.latex", preamble = "\\usepackage{xcolor}")

ax = fig.add_axes([0.1,0.24,0.89,0.74])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ylabel("\$[p.u.]\$")
xlabel("\$u\$")

colors = ["lightseagreen", "seagreen", "navy", "blue"]
#labls = ["\$\\hat{\\chi}\$", "\$\\chi\$", "\$\\chi^{+}\$", "\$\\chi^{-}\$"]
#ylim(bottom = 0.0, top = 80.0)
#ax.get_yaxis().set_visible(false)

x = [i for i in 1:length(uRESs)]
labls = [string(i) for i in 1:length(uRESs)]
xticks(x, labls)

plot(x, μ, color = colors[1], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\mu\$")
plot(x, μm, color = colors[2], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\mu^{+}\$")
plot(x, μp, color = colors[3], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\mu^{-}\$")

# for (i, val) in enumerate([χ0, χ1, χm, χp])
#     annotate(string(round(val, digits = 1)), (i - 0.4, val + 2.0))
# end

legend(loc = "lower right", fancybox = false, edgecolor = "black", framealpha = 0.85, ncol = 3)
savefig(string("plots_final//params_mu.pdf"), format = :pdf)

#-----------------
#-----------------

fig = figure(figsize=(5.0, 1.5))
rc("font", family = "serif", style = "italic", size = 11)
rc("text", usetex = true)
rc("lines", linewidth = 1)
rc("text.latex", preamble = "\\usepackage{xcolor}")

ax = fig.add_axes([0.1,0.24,0.89,0.74])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ylabel("\$[p.u.]\$")
xlabel("\$u\$")

colors = ["lightseagreen", "seagreen", "navy", "blue"]
#labls = ["\$\\hat{\\chi}\$", "\$\\chi\$", "\$\\chi^{+}\$", "\$\\chi^{-}\$"]
#ylim(bottom = 0.0, top = 80.0)
#ax.get_yaxis().set_visible(false)

x = [i for i in 1:length(uRESs)]
labls = [string(i) for i in 1:length(uRESs)]
xticks(x, labls)

plot(x, σ, color = colors[1], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\sigma\$")
plot(x, σm, color = "navy", marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\sigma^{+ \\backslash -}\$")

# for (i, val) in enumerate([χ0, χ1, χm, χp])
#     annotate(string(round(val, digits = 1)), (i - 0.4, val + 2.0))
# end

legend(loc = "upper left", fancybox = false, edgecolor = "black", framealpha = 0.85, ncol = 2)
savefig(string("plots_final//params_sigma.pdf"), format = :pdf)

#-----------------
#-----------------

fig = figure(figsize=(5.0, 1.5))
rc("font", family = "serif", style = "italic", size = 11)
rc("text", usetex = true)
rc("lines", linewidth = 1)
rc("text.latex", preamble = "\\usepackage{xcolor}")

ax = fig.add_axes([0.1,0.24,0.89,0.74])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ylabel("\$[p.u.]\$")
xlabel("\$u\$")

colors = ["lightseagreen", "seagreen", "navy", "blue"]
#labls = ["\$\\hat{\\chi}\$", "\$\\chi\$", "\$\\chi^{+}\$", "\$\\chi^{-}\$"]
#ylim(bottom = 0.0, top = 80.0)
#ax.get_yaxis().set_visible(false)

x = [i for i in 1:length(uRESs)]
labls = [string(i) for i in 1:length(uRESs)]
xticks(x, labls)

plot(x, σ, color = colors[1], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\sigma\$")
plot(x, σm, color = "navy", marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\sigma^{+ \\backslash -}\$")

# for (i, val) in enumerate([χ0, χ1, χm, χp])
#     annotate(string(round(val, digits = 1)), (i - 0.4, val + 2.0))
# end

legend(loc = "upper left", fancybox = false, edgecolor = "black", framealpha = 0.85, ncol = 2)
savefig(string("plots_final//params_sigma.pdf"), format = :pdf)

#-----------------
#-----------------

fig = figure(figsize=(5.0, 1.5))
rc("font", family = "serif", style = "italic", size = 11)
rc("text", usetex = true)
rc("lines", linewidth = 1)
rc("text.latex", preamble = "\\usepackage{xcolor}")

ax = fig.add_axes([0.1,0.24,0.89,0.74])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#ylabel("\$[p.u.]\$")
xlabel("\$u\$")

colors = ["lightseagreen", "seagreen", "navy", "blue"]
#labls = ["\$\\hat{\\chi}\$", "\$\\chi\$", "\$\\chi^{+}\$", "\$\\chi^{-}\$"]
ylim(bottom = -0.1, top = 0.8)
#ax.get_yaxis().set_visible(false)

x = [i for i in 1:length(uRESs)]
labls = [string(i) for i in 1:length(uRESs)]
xticks(x, labls)

#plot(x, χ0u ./ χ0, color = colors[1], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\beta_{u}\$")
plot(x, χ1u ./ χ1, color = colors[2], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\beta_{u}\$")
plot(x, χ1up ./ χ1, color = colors[3], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\beta^{+}_{u}\$")
plot(x, χ1um ./ χ1, color = colors[4], marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\beta^{-}_{u}\$")
#plot(x, χ1up ./ χ1, color = "navy", marker = "D", linewidth = 0.8, linestyle = "dashed", ms = 2.8, label = "\$\\sigma^{+ \\backslash -}\$")

# for (i, val) in enumerate([χ0, χ1, χm, χp])
#     annotate(string(round(val, digits = 1)), (i - 0.4, val + 2.0))
# end

legend(loc = "upper left", fancybox = false, edgecolor = "black", framealpha = 0.85, ncol = 3)
savefig(string("plots_final//beta.pdf"), format = :pdf)

χ1u
μ
