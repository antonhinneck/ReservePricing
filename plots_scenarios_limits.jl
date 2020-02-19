c10 = [5 / 255, 160 / 255, 250 / 255]
c20 = [5 / 255, 255 / 255, 50 / 255]
c1 = RGBA(c10..., 1.0)
c2 = RGBA(c20..., 1.0)
grad, grad_rgba, grad_rgb = gradient(c1, c2, length(limits))
styles = ["solid", "dashed", "dotted", "solid"]

## CHI MINUS
##----------

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.1,0.19,0.896,0.795])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$u\$")
ylabel("\$\\chi^{-}_{u}\$")
#xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

i = length(limits)
j = 1
for l in limits

    plot(string(u_buses), scenarios_chiM[j], color = grad_rgba[i], mec = grad_rgba[i], mfc = "white", label = "\$\\chi^{-}_{u}\$, \$$(l)\\overline{P}\$", lw = 1, ls = styles[j], marker = "D", ms = 3, mew = 1)
    global i -= 1
    global j += 1

end

legend(loc = "upper left", fancybox = false, edgecolor = "black")
savefig(string("plots_scenarios_limits//chi_m.pdf"), format = :pdf)

## CHI PLUS
##----------

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.1,0.19,0.896,0.795])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$u\$")
ylabel("\$\\chi^{+}_{u}\$")
#xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

i = length(limits)
j = 1
for l in limits

    plot(string(u_buses), -scenarios_chiP[j], color = grad_rgba[i], mec = grad_rgba[i], mfc = "white", label = "\$\\chi^{+}_{u}\$, \$$(l)\\overline{P}\$", lw = 1, ls = styles[j], marker = "D", ms = 3, mew = 1)
    global i -= 1
    global j += 1

end

legend(loc = "upper left", fancybox = false, edgecolor = "black")
savefig(string("plots_scenarios_limits//chi_p.pdf"), format = :pdf)

## DELTA P
##--------

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.1,0.19,0.896,0.795])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$i\$")
ylabel("\$\\overline{\\delta}_{i}\$")
#xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

styles = ["solid", "dashed", "dotted", "solid"]

i = length(limits)
j = 1
for l in limits

    plot([i for i in 1:n_generators], -scenarios_dP[j], color = grad_rgba[i], mec = grad_rgba[i], mfc = "white", label = "\$\\overline{\\delta}_{i}\$, \$$(l)\\overline{P}\$", lw = 1, ls = styles[j], marker = "D", ms = 3, mew = 1)
    global i -= 1
    global j += 1

end

legend(loc = "center right", fancybox = false, edgecolor = "black")
savefig(string("plots_scenarios_limits//delta_p.pdf"), format = :pdf)

## DELTA M
##--------

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.1,0.19,0.896,0.795])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$i\$")
ylabel("\$\\underline{\\delta}_{i}\$")
#xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

styles = ["solid", "dashed", "dotted", "solid"]

i = length(limits)
j = 1
for l in limits

    plot([i for i in 1:n_generators], -scenarios_dM[j], color = grad_rgba[i], mec = grad_rgba[i], mfc = "white", label = "\$\\underline{\\delta}_{i}\$, \$$(l)\\overline{P}\$", lw = 1, ls = styles[j], marker = "D", ms = 3, mew = 1)
    global i -= 1
    global j += 1

end

legend(loc = "center right", fancybox = false, edgecolor = "black")
savefig(string("plots_scenarios_limits//delta_m.pdf"), format = :pdf)
