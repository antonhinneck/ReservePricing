using Distributions
using PyPlot

sigma = 0.07
ϵ = 0.01
sn = Normal()
n = Normal(0, sigma)
tn = truncated(n, 0, Inf)

z0 = quantile(sn, 1 - ϵ)
z0 * 0.07
z1 = quantile(n, 1 - ϵ)
z2 = quantile(tn, 1 - 2 * ϵ)

sqrt(var(tn))
sqrt(var(n))

x = [i for i in range(-0.25,stop=0.25, length=80)]
yn = [pdf(n, i) for i in x]
ytn = [pdf(tn, i) for i in x]

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
xlabel("\$\\gamma\$")
ylabel("\$z_{u_{\\gamma\\sigma}}\$")

ax.axvline(z1)
plot(x, yn, color = "blue", mec = "white", mfc = "white", lw = 2, ls = "solid")
plot(x, ytn, color = "black", mec = "white", mfc = "white", lw = 1, ls = "dotted")
#plot(x, scenarios_zu, color = "black", mec = "white", mfc = "white", lw = 0.5, ls = "dashed")

#legend(loc = "upper left", fancybox = false, edgecolor = "black")
savefig(string("distribution.pdf"), format = :pdf)
