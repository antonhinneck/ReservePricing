fig = figure(figsize=(7, 2.0))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.128,0.265,0.868,0.72])
grid(linewidth = 0.2, linestyle = (0, (12, 4)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ylabel("Error \$[\\frac{\\\$}{h}],~[\\%]\$")
xlabel("Linear segments \$[\\frac{1}{MW p.u.}]\$")

yscale("log")

# quad_x = [range(coords[1][1], coords[length(coords)][1], step = 0.05);]
# quad_y = [quadraticCosts(generators, val, 5) for val in quad_x]
#
# c10 = [5 / 255, 160 / 255, 250 / 255]
# c20 = [5 / 255, 130 / 255, 200 / 255]
# c1 = RGBA(c10..., 1.0)
# c2 = RGBA(c20..., 1.0)
# grad, grad_rgba, grad_rgb = gradient(c1, c2, length(coords))
#
# for i in 1:length(coefs)
#     if true
#         x1 = coords[1][1]
#         x2 = coords[length(coords)][1]
#         y1 = coefs[i][1] * x1 + coefs[i][2]
#         y2 = coefs[i][1] * x2 + coefs[i][2]
#         if i == 1
#             plot([x1, x2],[y1, y2], color = grad_rgba[i], linewidth = 0.2, label = "Linear Approximation")
#         else
#             plot([x1, x2],[y1, y2], color = grad_rgba[i], linewidth = 0.2)
#         end
#     end
# end

plot([string(i) for i in 1:16], results_approx .- z1, ms = 1.8, marker = "D", color = "navy", linestyle = "dotted", linewidth = 0.5, label = "Error abs.")
plot([string(i) for i in 1:16], (results_approx ./ z1 .- 1.0) .* 100.0, ms = 1.8, marker = "D", color = "blue", linestyle = "dotted", linewidth = 0.5, label = "Error rel.")
legend(loc = "lower left", fancybox = false, edgecolor = "black")
savefig(string("error.pdf"), format = :pdf)
