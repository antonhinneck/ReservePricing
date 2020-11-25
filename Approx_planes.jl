function meshgrid(a::Array{T, 1} where T <: Real, b::Array{T, 1} where T <: Real)

	a_mesh = Array{Array{Float64, 1}, 1}(undef, length(b))
	for i in 1:length(b)
		#a_mesh[i] = Array{Float64, 1}(undef, length(a))
		a_mesh[i] = a
	end

	b_mesh = Array{Array{Float64}, 1}(undef, length(b))
	for i in 1:length(b)
		b_mesh[i] = Array{Float64, 1}(undef, length(a))
		for j in 1:length(a)
			b_mesh[i][j] = b[i]
		end
	end

	return a_mesh, b_mesh
end

function mesh_surface(a, b, f)

	z = Array{Array{Float64, 1}, 1}(undef, size(a, 1))

	for i in 1:size(a, 1)
		z[i] = Array{Float64, 1}(undef, size(a[1], 1))
		for j in 1:size(a[1], 1)
			z[i][j] = f(a[i][j], b[i][j])
		end
	end

	return z
end

function to2D(a)

	mat = Array{Float64, 2}(undef, (size(a, 1), size(a[1], 1)))

	for i in 1:size(a, 1)
		mat[i, :] = a[i]
	end

	return mat
end

# using PyPlot
# using Colors

#-------------------------------------
#-------------------------------------

c_vec_1 = [g.pi1  for g in generators]
c_vec_2 = [g.pi2  for g in generators]

p_min = [g.Pgmin  for g in generators]
p_max = [g.Pgmax  for g in generators]

g = 5

α_min = 0.0
α_max = 1.0

h1 = (x,y) -> @. c_vec_2[g] * (x * y) * (2.5)
h2 = (x,y) -> @.c_vec_2[g] * (x * y) * (-2.5)

env1 = (α,p) -> @.c_vec_2[g] * (α * p_max[g] + p * α_max - α_max * p_max[g]) * (2.5)
env2 = (α,p) -> @.c_vec_2[g] * (α_min * p + α * p_min[g] - α_min * p_min[g] ) * (2.5)
env3 = (α,p) -> @.c_vec_2[g] * (α_min * p + α * p_max[g] - α_min * p_max[g] ) * (2.5)
env4 = (α,p) -> @.c_vec_2[g] * (α_max * p + α * p_min[g] - α_max * p_min[g] ) * (2.5)

α = 0:0.1:1
p = 0:0.1:p_max[g]

using Gaston

surf(α, p, h2, w = :pm3d,
     Axes(title    = :Cost_Functions,
          palette  = :acton,
          #hidden3d = :on,
		  color = :blue,
		  view = "60, 30, 1.0, 1.0"))

surf!(α, p, h1, w = :pm3d)

surf(α, p, env2, w = :pm3d,
     Axes(title    = :Cost_Functions,
          palette  = :acton,
          #hidden3d = :on,
		  color = :blue,
		  view = "60, 230, 1.0, 1.0"))

surf!(α, p, h1, w = :pm3d)
surf!(α, p, env1, w = :pm3d)
surf!(α, p, env3, w = :pm3d)
surf!(α, p, env4, w = :pm3d)

# αv = collect(range(α_min; stop = α_max, step = 0.01))
# p = collect(range(p_min[g]; stop = p_max[g], step = 0.01))
# α = zeros(length(p))
# α[1:length(αv)] = αv

# f1(x,y) = c_vec_1[g] * (x^2 + y^2) + c_vec_2[g] * (x * y) * (1.0)
# f2(x,y) = c_vec_1[g] * (x^2 + y^2) * (sqrt(c_vec_2[g])/8)
# f3(x,y) = c_vec_1[g] * (x^2 + y^2)

g1 = (x,y) -> @. c_vec_1[g] * (x^2 + y^2)

g2 = (x,y) -> @. c_vec_1[g] * (x^2 + y^2) + c_vec_2[g] * (x * y) * (2.5)

g3 = (x,y) -> @. c_vec_1[g] * (x^2 + y^2) + c_vec_2[g] * (α_max * p_max[g]) * (2.5)

g4 = (x,y) -> @. c_vec_1[g] * (x^2 + y^2) + c_vec_2[g] * (x * y) * -(2.5)

g5 = (x,y) -> @. c_vec_1[g] * (x^2 + y^2) + c_vec_2[g] * (α_max * p_max[g]) * -(2.5)

#-------------------------------------
#-------------------------------------

surf(α, p, g5, w = :pm3d,
     Axes(title    = :Cost_Functions,
          palette  = :acton,
          #hidden3d = :on,
		  color = :blue,
		  view = "60, 230, 1.0, 1.0"))

surf!(α, p, g4, w = :pm3d)

surf!(α, p, g1, w = :pm3d)

surf!(α, p, g2, w = :pm3d)

surf!(α, p, g3, w = :pm3d)

#----------------
#----------------

surf(α, p, g5, w = :pm3d,
     Axes(title    = :Cost_Functions,
          palette  = :acton,
          #hidden3d = :on,
		  color = :blue,
		  view = "60, 230, 1.0, 1.0"))

surf!(α, p, g4, w = :pm3d)

surf!(α, p, g1, w = :pm3d)

surf!(α, p, g2, w = :pm3d)

surf!(α, p, g3, w = :pm3d)

#----------------
#----------------

surf(α, p, g4, w = :pm3d,
     Axes(title    = :Cost_Functions,
          palette  = :acton,
          hidden3d = :on,
		  view = "50, 30, 1.0, 1.0"))

surf!(α, p, g1, w = :pm3d)

surf!(α, p, g2, w = :pm3d)


#-------------------------------------
#-------------------------------------

fig = figure(figsize=(8, 8))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)
rc("axes", labelpad = 15)
#PyPlot.Axes3D()
#rc("axes", margins = 15)
#rc("axes.ticks.major", pad = 10)
# rc("ytick.major", pad = 10)
# rc("ztick.major", pad = 10)
# rc(ytick.major.pad = 80)
# rc(ztick.major.pad = 80

ylim(bottom = 0, top = 10000)
fig.add_axes(projection = "3d")

αv = collect(range(α_min; stop = α_max, step = 0.01))
p = collect(range(p_min[g]; stop = p_max[g], step = 0.01))
α = zeros(length(p))
α[1:length(αv)] = αv

f1(x,y) = c_vec_1[g] * (x^2 + y^2) + c_vec_2[g] * (x * y) * (1.0)
f2(x,y) = c_vec_1[g] * (x^2 + y^2) * (sqrt(c_vec_2[g])/8)
f3(x,y) = c_vec_1[g] * (x^2 + y^2)

zlabel("\$c_{i}\$")
ylabel("\$\\alpha_{u}\$")
xlabel("\$p_{i}\$")

a, b = meshgrid(α, p)
z1 = mesh_surface(a, b, f1)
z2 = mesh_surface(a, b, f2)
z3 = mesh_surface(a, b, f3)

plot_surface(to2D(b), to2D(a), to2D(z1), color = "orange", antialiased = false, rstride = 2, cstride = 2)
plot_surface(to2D(b), to2D(a), to2D(z2), color = "blue", antialiased = false, rstride = 2, cstride = 2)
plot_surface(to2D(b), to2D(a), to2D(z3), color = "lightgreen", antialiased = false, rstride = 2, cstride = 2)

gca()[:view_init](10, 140)
savefig(string("plots_final//test.png"), format = :png)

#------------------------------------------
#------------------------------------------

function Plot3D_costFunction(a, b, angle, id; dir = "mixed")

	fig = figure(figsize=(8, 8))
	rc("font", family = "serif", style = "italic", size = 14)
	rc("text", usetex = true)
	rc("lines", linewidth = 1)
	fig.add_axes(projection = "3d")

	f1(x,y) = c_vec_1[g] * (x^2 + y^2) #+ c_vec_2[g] * (x * y) * (1.0)
	f2(x,y) = c_vec_1[g] * (x^2 + y^2) + c_vec_2[g] * (x * y) * (1.0)
	f3(x,y) = c_vec_1[g] * (x^2 + y^2) + c_vec_2[g] * (x * y) * (-1.0)

	zlabel("\$c_{i}\$")
	ylabel("\$\\alpha_{u}\$")
	xlabel("\$p_{i}\$")

	z1 = mesh_surface(a, b, f1)
	z2 = mesh_surface(a, b, f2)
	z3 = mesh_surface(a, b, f3)

	plot_surface(to2D(b), to2D(a), to2D(z1), color = "orange", antialiased = true, rstride = 5, cstride = 5)
	plot_surface(to2D(b), to2D(a), to2D(z2), color = "green", antialiased = true, rstride = 5, cstride = 5)
	plot_surface(to2D(b), to2D(a), to2D(z3), color = "blue", antialiased = true, rstride = 5, cstride = 5)

    gca()[:view_init](15, angle)
    savefig(string("plots_final//$(dir)//3d_func$(id).png"), format = :png)

end

angles = [9.4 * i for i in 1:40]
for (i, angle) in enumerate(angles)
    Plot3D_costFunction(a, b, angle, i)
end

#--------------------------------------
#--------------------------------------

1 = 1

# using PyCall
# PyCall.current_python()

#ax = axes(projection='3d')
# z1 = c_vec_1[g] .* (α.^2 .+ p.^2) + c_vec_2[g] .* (p .* α)
# z2 = c_vec_1[g] .* (α.^2 .+ p.^2) + c_vec_2[g] .* (p .* α) .* (4.0)
#z = c_vec_1[g] ..* (a..^2 + b..^2) ..+ c_vec_2[g] ..* (a ..* b)

# Plots.plot(α, p, f, st=:surface, camera=(35, 10))
# Plots.savefig(string("plots_final//3d_func.pdf"))
# for (i, angl) in enumerate([9.5 * i for i in 1:40])
# 	Plots.plot(α, p, f, st=:surface, camera=(angl, 15))
# 	Plots.savefig(string("plots_final//positive//3d_func_$(i).png"))
# end


# z = Array{Float64, 2}(undef, (length(p), length(p)))
# #
# for (i, val1) in enumerate(α)
# 	for (j, val2) in enumerate(p)
# 		z[i, j] = f(val1, val2)
# 	end
# end

# #ax = fig.add_axes([0.095,0.19,0.90,0.795])
# grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
# ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

# ax.set_yscale("log")
# ax.set_axisbelow(true)

# xlim(left = 1, right = n_generators)
# zlim(bottom = -1.00, top = 3.5)
# ylim(bottom = 0.00, top = 1.5)
# xlim(left = 0.00, right = 1.0)

# plot_surface(a,b,z)
# plot3D((α,p),z,st=:surface,camera=(-30,30))
# surf(α, p, z1)
# surf(α, p, z2)

# legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9)

# m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG=output_level))
#
# ## Variables
# ##----------
# @variable(m, p[1:n_generators] >= 0)
# @variable(m, f[1:n_lines])
# @variable(m, θ[1:n_buses])
#
# ## General Constraints
# ##--------------------
# @constraint(m, θ[slack_bus] == 0)
# @expression(m, p_by_bus[i=1:n_buses], length(buses[i].genids) > 0 ? sum(p[k] for k in buses[i].genids) : 0.0)
# @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farmids) > 0 ? sum(farms[k].μ for k in buses[i].farmids) : 0.0)
# @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)
#
# @constraint(m, B * θ .== f)
# @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].u)
# @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].u)
#
# @constraint(m, cc1[i in 1:n_generators], p[i] <= generators[i].Pgmax)
# @constraint(m, cc2[i in 1:n_generators], -p[i] <= -generators[i].Pgmin)
#
# ## Generation Cost
# ##----------------
#
# @variable(m, d_con >= 0)
# @variable(m, d_lin >= 0)
# @variable(m, d_quad >= 0)
# @constraint(m, d_con == sum(generators[i].pi3 for i in 1:n_generators))
# @constraint(m, d_lin == sum(p[i] * generators[i].pi2 for i in 1:n_generators))
# @constraint(m, vec(vcat(0.5, d_quad, C_rt * p)) in RotatedSecondOrderCone())
# @expression(m, det_c, d_con + d_lin + d_quad)
#
# ## Objective
# ##----------
# @objective(m, Min, det_c)
#
# using PyPlot
#
# ####################
# ##  Prepare Data  ##
# ####################
# u = range(0.0,stop=2pi,length=300);
# v = range(0.0,stop=pi,length=300);
#
# lu = length(u);
# lv = length(v);
#
# x = zeros(lu,lv);
# y = zeros(lu,lv);
# z = zeros(lu,lv);
#
# for uu=1:lu
# 	for vv=1:lv
# 		x[uu,vv]= cos(u[uu])*sin(v[vv]);
# 		y[uu,vv]= sin(u[uu])*sin(v[vv]);
# 		z[uu,vv]= cos(v[vv]);
# 	end
# end
#
# #######################
# ##  Generate Colors  ##
# #######################
# colors = rand(lu,lv,3)
#
# ############
# ##  Plot  ##
# ############
# surf(x,y,z,facecolors=colors);
# savefig(string("plots_final//test.png"), format = :png)
