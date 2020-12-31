using Distributions, PyPlot

cd(@__DIR__)

mutable struct fe
    mu::Float64
    var::Float64
end

fes =  [[-0.0036, 0.0343],
        [0.0012, 0.0389],
        [-0.0141, 0.0381],
        [0.0002, 0.0405],
        [-0.0023, 0.0438],
        [0.0122, 0.0555],
        [0.0002, 0.0405],
        [0.0163, 0.0452],
        [-0.0, 0.0007],
        [-0.0003, 0.0065],
        [-0.0021, 0.0128],
        [-0.0006, 0.0037],
        [0.0109, 0.0443],
        [-0.0029, 0.0258]]

a = [i for i in 1:length(fes)]

caps = [70.0, 147.0, 102.0, 105.0, 113.0, 84.0, 59.0, 250.0, 118.0, 76.0, 72.0]

dn = Distributions.Normal()

x = [i for i in -6:0.01:6]
y = [pdf(dn, i) for i in x]

fes_farms = Vector{Vector{Float64}}()

for (i, farm) in enumerate(farms)
        push!(fes_farms, [fes[i][1] * caps[i] / 100, fes[i][2]^2 * (caps[i] / 100)^2])
end

using Colors

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

scenarios = [string(i) for i in 1:length(fes_farms)]
mg = gradient(RGBA(15 / 255, 25 / 255, 120 / 255, 1.0), RGBA(0 / 255, 255 / 255, 80 / 255, 1.0), length(fes_farms))

g_rgba = Vector{Tuple{Float64, Float64, Float64, Float64}}()
for i in 1:length(mg)
    push!(g_rgba, ((mg[i].r,mg[i].g,mg[i].b,1.0)))
end

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.095,0.19,0.9,0.795])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

xlabel("\$u\$")
ylabel("\$\\chi^{\\pm}_{(u)},~\\beta^{\\pm}\$")
ylim(bottom = 0.0, top = 0.15)

c = 1
for farm in fes_farms
        myd = Normal(farm[1], sqrt(farm[2]))
        x = [i for i in -6:0.01:6]
        y = [pdf(myd, i) for i in x]
        plot(x, y, color = g_rgba[c], lw = 0.5, ls = "solid")
        global c += 1
end

plot(x, y, color = "black", lw = 0.8, ls = "solid")

#legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9)
savefig(string("dist.pdf"), format = :pdf)

for farm in fes_farms
        print(string(farm[1], " & "))
end
println()
for farm in fes_farms
        print(string(sqrt(farm[2]), " & "))
end
