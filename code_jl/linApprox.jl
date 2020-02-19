function quadraticCosts(generators::Vector{Generator}, x::Float64, g::Int64)

    return generators[g].pi1 * x * x + generators[g].pi2 * x + generators[g].pi3

end

mutable struct aprx
    coords::Vector{Tuple{Float64, Float64}}
    coefs::Vector{Tuple{Float64, Float64}}
end

function approx(generators::Vector{Generator}, g::Int64; error = 10e-9)

    left = generators[g].Pgmin
    right = generators[g].Pgmax
    segments = ceil(right) * 10
    stepsize = (right - left) / segments
    coords = Vector{Tuple{Float64, Float64}}()
    coefs = Vector{Tuple{Float64, Float64}}()
    for i in 1:segments
        x0 = left + stepsize * (i - 1)
        x1 = left + stepsize * i
        y0 = quadraticCosts(generators, x0, g)
        y1 = quadraticCosts(generators, x1, g)
        m = (y1 - y0) / (x1 - x0)
        n = y1 - m * x1
        @assert abs(n - (y0 - m * x0)) <= error
        push!(coefs, (m,n))
        if i == 1
            push!(coords, (x0,x0))
        end
        push!(coords, (x1,y1))
    end
    return coords, coefs
end

my_aprxs = Vector{aprx}()
for i in 1:n_generators
    push!(my_aprxs, aprx(approx(generators, i)...))
end
my_aprxs
coords = my_aprxs[5].coords
coefs = my_aprxs[5].coefs
n_coefs = length(coefs)

fig = figure(figsize=(6, 5))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.15,0.09,0.84,0.9])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ylabel("\$c(p_{5})\$")
xlabel("\$p_{5}\$")
#xlim(left = 0, right = 0.6)
#ylim(bottom = -100, top = 600)

quad_x = [range(coords[1][1], coords[length(coords)][1], step = 0.05);]
quad_y = [quadraticCosts(generators, val, 5) for val in quad_x]

c10 = [5 / 255, 160 / 255, 250 / 255]
c20 = [5 / 255, 160 / 255, 250 / 255]
#c20 = [0 / 255, 100 / 255, 255 / 255]
c1 = RGBA(c10..., 1.0)
c2 = RGBA(c20..., 1.0)
grad, grad_rgba, grad_rgb = gradient(c1, c2, length(coords))
#cm = PyPlot.LinearSegmentedColormap.from_list("Linear Approximations", [Tuple(c10), Tuple(c20)], length(grad_rgba))

xdescent = 0
ydescent = 0
width = 4
height = 1
num_stripes = length(grad_rgba)
stripes = Vector{Any}()
for i in 1:num_stripes
    s = PyPlot.matplotlib.patches.Rectangle([xdescent + i * width / num_stripes, ydescent],
                          width / num_stripes,
                          height,
                          fc=cm((2 * i + 1) / (2 * num_stripes)))
            push!(stripes, s)
end

for i in 1:length(coefs)
    x1 = coords[1][1]
    x2 = coords[length(coords)][1]
    y1 = coefs[i][1] * x1 + coefs[i][2]
    y2 = coefs[i][1] * x2 + coefs[i][2]
    if i == 1
        plot([x1, x2],[y1, y2], color = grad_rgba[i], linewidth = 0.2, label = "Linear Approximation")
    else
        plot([x1, x2],[y1, y2], color = grad_rgba[i], linewidth = 0.2)
    end
end

plot(quad_x, quad_y, color = "black", linestyle = "solid", linewidth = 1, label = "Quadratic Costs")
legend(loc = "upper left", fancybox = false, edgecolor = "black")
savefig(string("approx.pdf"), format = :pdf)
