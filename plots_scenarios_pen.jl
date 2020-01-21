function gradient(start::RGBA, stop::RGBA, n_el::T where T <: Integer)

    grad = Vector{RGB}()

    stepsize_r = (stop.r - start.r) / n_el
    stepsize_g = (stop.g - start.g) / n_el
    stepsize_b = (stop.b - start.b) / n_el

    push!(grad, start)

    current = [start.r, start.g, start.b, 1.0]
    for i in 1:(n_el - 1)
        current[1] += stepsize_r
        current[2] += stepsize_g
        current[3] += stepsize_b
        push!(grad, RGBA(current...))
    end

    return grad
end

mg = gradient(RGBA(15 / 255, 25 / 255, 120 / 255, 1.0), RGBA(0 / 255, 255 / 255, 80 / 255, 1.0), length(scalings))

g_rgba = Vector{Tuple{Float64, Float64, Float64, Float64}}()
for i in 1:length(mg)
    push!(g_rgba, ((mg[i].r,mg[i].g,mg[i].b,1.0)))
end

fig = figure(figsize=(8, 3.6))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.0925,0.12,0.9,0.86])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$u\$")
ylabel("\$\\chi^{(-)}_{u}\$")
##xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

for i in 1:length(mg)

    j = length(mg) - (i - 1)
    #plot(u_buses, scenarios_sigma[j], color = [mg[j].r,mg[j].g,mg[j].b, 1.0], mec = [mg[j].r,mg[j].g,mg[j].b, 1.0], mfc = "white", lw = 1, ls = "dashed", marker = "D", ms = 3, mew = 1) #label = "\$\\chi_{u}\$, \$$(scalings[i])\\sigma_{u}\$",
    plot(u_buses, scenarios_chi[i], color = [mg[j].r,mg[j].g,mg[j].b, 1.0], mec = [mg[j].r,mg[j].g,mg[j].b, 1.0], mfc = "white", label = "\$$(scalings[i])\\sigma_{u}\$, \$\\chi^{-}_u\$", lw = 1, ls = "dotted", marker = "D", ms = 3, mew = 1)

end

legend(loc = "upper left", fancybox = false, edgecolor = "black")
savefig(string("plots_scenarios_pen//s_chi.pdf"), format = :pdf)

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
##xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

plot(scenarios, scenarios_zu, color = "black", mec = "white", mfc = "white", lw = 0.5, ls = "dashed") #label = "\$\\chi_{u}\$, \$$(scalings[i])\\sigma_{u}\$",
for i in 1:length(mg)
    plot(scenarios[i], scenarios_zu[i], color = "white", mec = g_rgba[i], mfc = "white", lw = 0, marker = "D", ms = 3.6, mew = 1.6)
end

#legend(loc = "upper left", fancybox = false, edgecolor = "black")
savefig(string("plots_scenarios_pen//s_zu.pdf"), format = :pdf)

fig = figure(figsize=(4, 3.6))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.185,0.13,0.8075,0.8625])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$\\gamma\$")
ylabel("\$z_{\\gamma\\sigma}\$")
##xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

plot(scenarios, scenarios_z, color = "black", mec = "white", mfc = "white", lw = 0.5, ls = "dashed") #label = "\$\\chi_{u}\$, \$$(scalings[i])\\sigma_{u}\$",
for i in 1:length(mg)
    plot(scenarios[i], scenarios_z[i], color = "white", mec = g_rgba[i], mfc = "white", lw = 0, marker = "D", ms = 3.6, mew = 1.6)
end

#legend(loc = "upper left", fancybox = false, edgecolor = "black")
savefig(string("plots_scenarios_pen//s_z.pdf"), format = :pdf)
