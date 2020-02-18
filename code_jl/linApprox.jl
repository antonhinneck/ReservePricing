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
    segments = ceil(right) * 8
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
coords = my_aprxs[1].coords
coefs = my_aprxs[1].coefs
n_coefs = length(coefs)

fig = figure(figsize=(8, 6))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.09,0.09,0.9,0.9])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ylabel("\$c(p_{1})\$")
xlabel("\$p_{1}\$")
xlim(left = 0, right = 0.6)
#ylim(bottom = -100, top = 600)

quad_x = [range(coords[1][1], 1, step = 0.01);]
quad_y = [quadraticCosts(generators, val, 1) for val in quad_x]

clrs = ["orange","yellow","lightgreen","green","teal", "blue", "orange","yellow","lightgreen","green","teal", "blue",
        "orange","yellow","lightgreen","green","teal", "blue", "orange","yellow","lightgreen","green","teal", "blue",
        "orange","yellow","lightgreen","green","teal", "blue", "orange","yellow","lightgreen","green","teal", "blue",
        "orange","yellow","lightgreen","green","teal", "blue", "orange","yellow","lightgreen","green","teal", "blue"]
for i in 1:length(coefs)
    x1 = coords[1][1]
    x2 = coords[length(coords)][1]
    y1 = coefs[i][1] * x1 + coefs[i][2]
    y2 = coefs[i][1] * x2 + coefs[i][2]
    #println(string(x1," ",y1," ",x2," ",y2))
    plot([x1, x2],[y1, y2], color = clrs[i], linewidth = 0.2, label = string("segment ",i))
end

plot(quad_x, quad_y, color = "black", linestyle = "dotted", linewidth = 1.5, label = "Quadratic Costs")

legend(loc = "upper left", fancybox = false, edgecolor = "black")
savefig(string("approx.pdf"), format = :pdf)
