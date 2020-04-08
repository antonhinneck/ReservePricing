dist_parent = Distributions.Normal(0,1)
dist_minus = Distributions.TruncatedNormal(dist_parent.μ, dist_parent.σ, 0, Inf64)
dist_plus = Distributions.TruncatedNormal(dist_parent.μ, dist_parent.σ,-Inf64, 0)

x = [i for i in range(-5.0,5.0, step = 0.01)]
xm = [i for i in range(-5.0,5.0, step = 0.01) if i >= 0]
xp = [i for i in range(-5.0,5.0, step = 0.01) if i < 0]
ω = [pdf(dist_parent, i) for i in x]
ωm = [pdf(dist_minus, i) for i in xm]
ωp = [pdf(dist_plus, i) for i in xp]

fig = figure(figsize=(4, 3.6))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.18,0.13,0.81,0.8625])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)
#ylim(bottom = 0, top = 2490)
#ax.set_yscale("log")
#ax.set_axisbelow(true)
ylabel("\$\\phi(\\omega^{(\\mp)})\$")
xlabel("\$\\omega\$")
##xlim(left=-5,right=5)

plot(x, ω, color = "navy", mec = "white", mfc = "white", lw = 1.0, ls = "solid", label = "\$\\omega\$") #label = "\$\\chi_{u}\$, \$$(scalings[i])\\sigma_{u}\$",
plot(xp, ωp, color = "teal", mec = "white", mfc = "white", lw = 1.0, ls = "dashed", label = "\$\\omega^{+}\$")
plot(xm, ωm, color = "lightgreen", mec = "white", mfc = "white", lw = 1.0, ls = "dashed", label = "\$\\omega^{-}\$")

legend(loc = "upper left", fancybox = false, edgecolor = "black")
savefig(string("plots_dists//truncated.pdf"), format = :pdf)
