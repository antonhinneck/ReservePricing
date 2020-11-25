cd(@__DIR__)
include("pkgs.jl")
include("code_jl/input.jl")
include("code_jl/gradients.jl")
include("code_jl/TruncatedGaussian.jl")

case_data = load("data//118bus.jld")
buses = case_data["buses"]
generators = case_data["generators"]

n_buses = size(buses, 1)
n_generators = size(generators, 1)

slack_bus = findall(b -> b.kind == :Ref, buses)[1]

lines = case_data["lines"]
line_limits = [ 175	175	500	175	175	175	500	500	500	175	175	175	175	175	175	175	175	175	175	175	500	175	175	175	175	175	175	175	175	175	500	500	500	175	175	500	175	500	175	175	140	175	175	175	175	175	175	175	175	500	500	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	500	175	175	500	500	500	500	500	500	500	175	175	500	175	500	175	175	500	500	175	175	175	175	175	175	175	500	175	175	175	175	175	175	500	500	175	500	500	200	200	175	175	175	500	500	175	175	500	500	500	175	500	500	175	175	175	175	175	175	175	175	175	175	200	175	175	175	175	175	175	175	175	175	500	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	500	175	175	175	500	175	175	175]
thermalLimitscale = 1
for i in 1:length(lines)
    lines[i].u = 1.0 * thermalLimitscale * line_limits[i] / 100 #0.99
end
n_lines = size(lines, 1)
x_vec = [l.x for l in lines]
X = diagm(0 => x_vec)

A = zeros(n_lines, n_buses)
for i in 1:n_buses
    if size(buses[i].outlist, 1) > 0
        for l in buses[i].outlist
            A[lines[l].arcID, i] = -1
        end
    end
    if size(buses[i].inlist, 1) > 0
        for l in buses[i].inlist
            A[lines[l].arcID, i] = 1
        end
    end
end

# for i in σ_vec
#     print(string(i,", "))
# end

B = X ^ (-1) * A
B_node = A' * B

d = [b.Pd for b in buses]
sum([b.Pd for b in buses])
sum([g.Pgmax for g in generators])

#include("code_jl/linApprox.jl")

include("code_jl//farms.jl")
farms, n_farms, σ_vec, Σ, s_sq, Σ_rt, s = create_wind_farms(scaling_sigma = 1.0)
σ_vec = [f.σ for f in farms]
u_buses = [f.bus for f in farms]
μ_vec = [f.μ for f in farms]
p_U = μ_vec
ν = sum(μ_vec)
sqrt(0.00178)
include("code_jl/TruncatedGaussian.jl")
lower, upper = splitGaussians(zeros(length(μ_vec)), [f.σ for f in farms], 0.0)
μm = upper[1]
Σm = upper[2]
Σm_rt = upper[3]
sm_sq = upper[4]
sm = sqrt(sm_sq)#sqrt(sm_sq)
μp = lower[1]
Σp = lower[2]
Σp_rt = lower[3]
sp_sq = lower[4]
sp = sqrt(sp_sq)
sum(Σp_rt)

muF = [μp..., -μp...]
ΣF = zeros(n_farms * 2,n_farms * 2)
ΣF[1:n_farms, 1:n_farms] = diagm(σ_vec)
ΣF[(n_farms + 1):(2 * n_farms), (n_farms + 1):(2 * n_farms)] = diagm(σ_vec)
ΣF

μm = [0 for i in farms]
μp = [0 for i in farms]

σm = Vector{Float64}()
for i in 1:size(Σm, 1)
    push!(σm, Σm[i,i])
end

ϵ = 0.01
d = Normal()
z = quantile(d, 1 - ϵ)
zCh = sqrt((1 - ϵ)/ϵ)
dT = TruncatedNormal(μm[1], Σm_rt[1,1], 0, Inf64)
Σm_rt[1,1] * z

zn = quantile(d, 1 - (ϵ / 0.5))

σ_vec
mysum = 0
for s in σ_vec
    dist = Normal(0,s)
    global mysum += quantile(dist, 1 - ϵ)
end

mysum
z
z * sum(σ_vec)
σ_vec = σ_vec
zt = Vector{Float64}()
μm = Vector{Float64}()
varm = Vector{Float64}()
mysum = 0
for i in 1:length(σ_vec)
    dist = truncated(Normal(0, σ_vec[i]),0,Inf64)
    #distp = truncated(Normal(0, σ_vec[i]),0,Inf64)
    global mysum += quantile(dist, 1 - (2 * ϵ))
    push!(zt, quantile(dist, 1 - (2 * ϵ)))
    push!(μm, mean(dist))
    push!(varm, var(dist))
end
μp = -μm
varp = varm

sum(μm)
n = Normal(2, 0)
nt = truncated(n,2,Inf)
mean(nt)

zz = quantile(truncated(Normal(),0,Inf), 1 - (2 * ϵ))

s
varm
mysum
sm
zz * sqrt(sum(Σm))
z * sqrt(varm)
sum(zt)
z
zz
Σm_rt

sum(μm)
z * sqrt(sum(Σm))
z * sqrt(sum(Σ))
zz * sqrt(sum(varm))

zm = diagm(zt)
sum(zm * Σm)
sum(zm) * sum(Σm)

#include("plotDists.jl")

d = [b.Pd for b in buses]

counter = 1
for f in farms
    print(string(counter," & "))
    global counter += 1
end

counter = 1
for f in farms
    print(string(round(f.σ, digits = 4)," & "))
    global counter += 1
end

counter = 1
for f in farms
    print(string(round(μm[counter], digits = 4)," & "))
    global counter += 1
end

counter = 1
for f in farms
    print(string(round(Σm_rt[counter,counter], digits = 4)," & "))
    global counter += 1
end

for (i,f) in enumerate(farms)
    push!(buses[f.bus].farmids, i)
end

## Stochastic parameters
########################

ϵ = 0.01
z = quantile(Normal(0,1), 1-ϵ)
za = quantile(Normal(0,1), 1 - ϵ / 0.5)

## Generation Costs
###################

c_vec = [g.pi1  for g in generators]
C_mat = diagm(0 => c_vec)
C_rt = sqrt(C_mat)

#include("Approx_planes.jl")

function updateGen(min::Float64, max::Float64)

    @assert min >= 0 && min <= 1 && max <= 1 && max >= 0 "Values must be between 0 and 1."

    case_data = load("data//118bus.jld")
    generators = case_data["generators"]

    for (i, g) in enumerate(generators)
         g.Pgmin = g.Pgmax * min
         g.Pgmax = g.Pgmax * max
    end
    return case_data, generators
end

μ = zeros(n_farms) .+ 0.1
μs = sum(μ)
𝛭 = sum(μ[f] for f in 1:n_farms)

ζ1 = (sqrt(pi / 2))
ζ2 = (2 * pi - 4) / (2 * pi)


μm = μ .+ (σ_vec .* (1 / ζ1))
μp = μ .- (σ_vec .* (1 / ζ1))
μms = sum(μm)
μps = sum(μp)

## Experiments
##############

## Linear Costs
###############
updateGen(0.2, 0.5)
# 0.1 0.39
#
# 105939 - 105157
#
#
# 0.1 0.395
#
# 104911 - 104878
#
# 0.1 0.4
#
# 104611 - 104603
#include("code_jl/linApprox.jl")

include("models/dccc.jl")
m_dccc = build_dccc(generators, buses, lines, farms)
optimize!(m_dccc)
termination_status(m_dccc)
z1 = objective_value(m_dccc)
value(m_dccc[:det_c])
value(m_dccc[:d_lin])
value(m_dccc[:d_con])
value(m_dccc[:d_quad])
value(m_dccc[:unc_c])

include("models/dccc_apx.jl")
mylim = 0.00001
m_dccc_apx = build_dccc_apx(generators, buses, lines, farms, α_min = [mylim for i in 1:n_generators])
optimize!(m_dccc_apx)
objective_value(m_dccc_apx)
sum(value.(m_dccc_apx[:p]))
value.(m_dccc_apx[:α])
sum(value.(m_dccc_apx[:α]))
sum(value.(m_dccc_apx[:ψ]))
value(m_dccc_apx[:det_c])
value(m_dccc_apx[:d_lin])
value(m_dccc_apx[:d_con])
value(m_dccc_apx[:d_quad])
value(m_dccc_apx[:unc_c])
value(m_dccc_apx[:d_bil])

include("models/dccc_det.jl")
m_dccc_det = build_dccc_det(generators, buses, lines, farms, value.(m_dccc_apx[:α]), value.(m_dccc_apx[:p]))
optimize!(m_dccc_det)
objective_value(m_dccc_det)
value(m_dccc_det[:d_lin])
value(m_dccc_det[:d_con])
value(m_dccc_det[:d_quad])
value(m_dccc_det[:unc_c])
value(m_dccc_det[:d_bil])

function adjust_limits_sw(alphas::Array{T, 1} where T <: Real; lim = Float64(0.0), dec = 6)

    new_min = zeros(length(alphas))
    new_max = zeros(length(alphas))

    for (i, a) in enumerate(alphas)
        new_min[i] = max(round(a - lim, digits = dec), 0.0)
        new_max[i] = min(round(a + lim, digits = dec), 1.0)
    end

    return new_min, new_max
end

function min_alg_sws(iters::Int64; lim_init = 0.00001, lim_const = 0.001, scheduler = 1.0, verbose = true)

    progression_apx = Vector{Float64}()
    progression_det = Vector{Float64}()
    iteration = Vector{String}()

    lims_min = nothing
    lims_max = nothing

    current_apx_model = nothing
    ctr = 0

    for i in 1:iters

        if i == 1
            m_dccc_apx_alg = build_dccc_apx(generators, buses, lines, farms, α_min = [lim_init for i in 1:n_generators])
            optimize!(m_dccc_apx_alg)
        else
            m_dccc_apx_alg = build_dccc_apx(generators, buses, lines, farms, α_min = lims_min, α_max = lims_max, output_level = 0)
            optimize!(m_dccc_apx_alg)
        end

        m_dccc_det_alg = build_dccc_det(generators, buses, lines, farms, value.(m_dccc_apx_alg[:α]), value.(m_dccc_apx_alg[:p]), output_level = 0)
        optimize!(m_dccc_det_alg)

        if verbose
            println(string("ITR: ", i, " APX: ", objective_value(m_dccc_apx_alg), " DET: ", objective_value(m_dccc_det_alg)))
        end

        lims_min, lims_max = adjust_limits_sw(value.(m_dccc_apx_alg[:α]), lim = lim_const)

        if i == 1
            push!(progression_apx, objective_value(m_dccc_apx_alg))
            push!(progression_det, objective_value(m_dccc_det_alg))
            push!(iteration, string(1))
            ctr += 1
        else
            if progression_det[ctr] > objective_value(m_dccc_det_alg) && termination_status(m_dccc_det_alg) == MOI.TerminationStatusCode(1) && objective_value(m_dccc_det_alg) != 0.0
                push!(iteration, string(i))
                push!(progression_apx, objective_value(m_dccc_apx_alg))
                push!(progression_det, objective_value(m_dccc_det_alg))
                current_apx_model = m_dccc_apx_alg
                ctr += 1
            end
        end

        #lim_const = lim_const * scheduler
    end

    return current_apx_model, [progression_apx, progression_det, iteration]
end

min_model, data = min_alg_sws(60, lim_init = 0.00, lim_const = 0.006)
objective_value(min_model)

m_dccc_det = build_dccc_det(generators, buses, lines, farms, value.(min_model[:α]), value.(min_model[:p]), output_level = 1)
optimize!(m_dccc_det)
objective_value(m_dccc_det)
value(m_dccc_det[:det_c])
value(m_dccc_det[:d_lin])
value(m_dccc_det[:d_con])
value(m_dccc_det[:d_quad])
value(m_dccc_det[:unc_c])
value(m_dccc_det[:d_bil])

using PyPlot

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.09,0.2,0.905,0.78])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("Iteration")
ylabel("\$z\$")
ylim(bottom = 98260, top = 98550)
#ylabel("\$\\chi^{+}_{u}\$")

plot(data[3], data[1], color = "lightgreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$z^{APX}\$")
plot(data[3], data[2], color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$z^{DET}\$")

legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9)
savefig(string("plots_final//algorithm.pdf"), format = :pdf)
#sm = sqrt(ζ2) * s

include("models/dccc_a_apx.jl")
m_dccc_a_apx = build_dccc_a_apx(generators, buses, lines, farms, αm_min = ones(n_generators) * 0.005, αp_min = ones(n_generators) * 0.005)
optimize!(m_dccc_a_apx)
termination_status(m_dccc_a_apx)
z4 = objective_value(m_dccc_a_apx)
value(m_dccc_a_apx[:det_c])
value(m_dccc_a_apx[:d_lin])
value(m_dccc_a_apx[:d_con])
value(m_dccc_a_apx[:d_quad])
value(m_dccc_a_apx[:unc_c])
value(m_dccc_a_apx[:d_bil])
value.(m_dccc_a_apx[:αm])
value.(m_dccc_a_apx[:αp])

include("models/dccc_a_det.jl")
m_dccc_a_det = build_dccc_a_det(generators, buses, lines, farms, value.(m_dccc_a_apx[:αm]), value.(m_dccc_a_apx[:αp]), value.(m_dccc_a_apx[:p]))
optimize!(m_dccc_a_det)
termination_status(m_dccc_a_det)
z4 = objective_value(m_dccc_a_det)
value(m_dccc_a_det[:det_c])
value(m_dccc_a_det[:d_lin])
value(m_dccc_a_det[:d_con])
value(m_dccc_a_det[:d_quad])
value(m_dccc_a_det[:unc_c])
value(m_dccc_a_det[:d_bil])

value.(m_dccc_det[:α])
value.(m_dccc_a_det[:αp])

function sw2n2n(alpha)

    # g x f
    new_alpha = ones(n_generators, n_farms)
    for i in 1:n_farms
        new_alpha[:, i] = [alpha[i] for i in 1:n_generators]
    end
    return new_alpha

end

α_min_init = sw2n2n(value.(m_dccc_a_det[:αp])) .- 0.0002 #ones((n_generators, n_farms)) * 0.000001
α_max_init = ones((n_generators, n_farms)) * 1.0
include("models/dccc_n2n_apx.jl")
m_dccc_n2n_apx = build_dccc_n2n_apx(generators, buses, lines, farms, α_min_init, α_max_init)
optimize!(m_dccc_n2n_apx)
termination_status(m_dccc_n2n_apx)
z4 = objective_value(m_dccc_n2n_apx)
value(m_dccc_n2n_apx[:det_c])
value(m_dccc_n2n_apx[:d_lin])
value(m_dccc_n2n_apx[:d_con])
value(m_dccc_n2n_apx[:d_quad])
value(m_dccc_n2n_apx[:unc_c])
value(m_dccc_n2n_apx[:u_bil])
value.(m_dccc_n2n_apx[:p_uncert])

α_det = value.(m_dccc_n2n_apx[:α])

include("models/dccc_n2n_det.jl")
m_dccc_n2n_det = build_dccc_n2n_det(generators, buses, lines, farms, α_det)
optimize!(m_dccc_n2n_det )
termination_status(m_dccc_n2n_det)
z4 = objective_value(m_dccc_n2n_det)
value(m_dccc_n2n_det[:det_c])
value(m_dccc_n2n_det[:d_lin])
value(m_dccc_n2n_det[:d_con])
value(m_dccc_n2n_det[:d_quad])
value(m_dccc_n2n_det[:unc_c])
value(m_dccc_n2n_det[:u_bil])

function adjust_limits_a(alphas::Array{T, 2} where T <: Real; lim = Float64(0.0), dec = 6)

    new_min = deepcopy(alphas)
    new_max = deepcopy(alphas)

    new_min -= ones(size(alphas)) * lim
    new_max += ones(size(alphas)) * lim

    return new_min, new_max
end

function min_alg_a(iters::Int64, α_initm, α_initp; lim_init = 0.01, lim_const = 0.01, scheduler = 1.0, verbose = true)

    progression_apx = Vector{Float64}()
    progression_det = Vector{Float64}()
    iteration = Vector{String}()

    # if α_initm == nothing
    #     α_initm = zeros(n_generators, n_farms)
    #     α_initp = zeros(n_generators, n_farms)
    # end


    lims_min_m = nothing
    lims_max_m = nothing

    lims_min_p = nothing
    lims_max_p = nothing

    current_apx_model = nothing
    ctr = 0

    for i in 1:iters

        if i == 1
            #m_dccc_n2n_a_apx_alg = build_dccc_n2n_a_apx(generators, buses, lines, farms, ones(n_generators, n_farms) * lim_init, ones(n_generators, n_farms) * 1.0, ones(n_generators, n_farms) * lim_init, ones(n_generators, n_farms) * 1.0
            m_dccc_n2n_a_apx_alg = build_dccc_n2n_a_apx(generators, buses, lines, farms, α_initm .- lim_init, ones(n_generators, n_farms), α_initp .- lim_init, ones(n_generators, n_farms))
            optimize!(m_dccc_n2n_a_apx_alg)
            println(objective_value(m_dccc_n2n_a_apx_alg))
        else
            m_dccc_n2n_a_apx_alg = build_dccc_n2n_a_apx(generators, buses, lines, farms, lims_min_m, lims_max_m, lims_min_p, lims_max_p)
            optimize!(m_dccc_n2n_a_apx_alg)
        end

        m_dccc_n2n_a_det_alg = build_dccc_n2n_a_det(generators, buses, lines, farms, value.(m_dccc_n2n_a_apx_alg[:αm]), value.(m_dccc_n2n_a_apx_alg[:αp]), output_level = 0)
        optimize!(m_dccc_n2n_a_det_alg)

        if verbose
            println(string("ITR: ", i, " APX: ", objective_value(m_dccc_n2n_a_apx_alg), " DET: ", objective_value(m_dccc_n2n_a_det_alg)))
        end

        lims_min_m, lims_max_m = adjust_limits_a(value.(m_dccc_n2n_a_apx_alg[:αm]), lim = lim_const)
        lims_min_p, lims_max_p = adjust_limits_a(value.(m_dccc_n2n_a_apx_alg[:αp]), lim = lim_const)

        if i == 1
            push!(progression_apx, objective_value(m_dccc_n2n_a_apx_alg))
            push!(progression_det, objective_value(m_dccc_n2n_a_apx_alg))
            push!(iteration, string(1))
            ctr += 1
        else
            if progression_det[ctr] > objective_value(m_dccc_n2n_a_det_alg) && termination_status(m_dccc_n2n_a_det_alg) == MOI.TerminationStatusCode(1) && objective_value(m_dccc_n2n_a_det_alg) != 0.0
                push!(iteration, string(i))
                push!(progression_apx, objective_value(m_dccc_n2n_a_apx_alg))
                push!(progression_det, objective_value(m_dccc_n2n_a_apx_alg))
                current_apx_model = m_dccc_n2n_a_apx_alg
                ctr += 1
            end
        end

        #lim_const = lim_const * scheduler
    end

    return current_apx_model, [progression_apx, progression_det, iteration]
end

# α_min_initm = ones((n_generators, n_farms)) * 0.0001
α_min_initm = sw2n2n(value.(m_dccc_a_det[:αm])) .- 0.001 # value.(m_dccc_n2n_det[:α]) # ones((n_generators, n_farms)) * 0.01
α_max_initm = ones((n_generators, n_farms)) * 1.0
# α_min_initp = ones((n_generators, n_farms)) * 0.0001
α_min_initp = sw2n2n(value.(m_dccc_a_det[:αp])) .- 0.001 # value.(m_dccc_n2n_det[:α]) # ones((n_generators, n_farms)) * 0.01
α_max_initp = ones((n_generators, n_farms)) * 1.0
include("models/dccc_n2n_a_apx.jl")
m_dccc_n2n_a_apx = build_dccc_n2n_a_apx(generators, buses, lines, farms, α_min_initm, α_max_initm, α_min_initp, α_max_initp)
optimize!(m_dccc_n2n_a_apx)
value(m_dccc_n2n_a_apx[:det_c])
value(m_dccc_n2n_a_apx[:d_lin])
value(m_dccc_n2n_a_apx[:d_con])
value(m_dccc_n2n_a_apx[:d_quad])
value(m_dccc_n2n_a_apx[:unc_c])
value(m_dccc_n2n_a_apx[:u_bil])

α_detm = value.(m_dccc_n2n_a_apx[:αm])
α_detp = value.(m_dccc_n2n_a_apx[:αp])

include("models/dccc_n2n_a_det.jl")
m_dccc_n2n_a_det = build_dccc_n2n_a_det(generators, buses, lines, farms, α_detm, α_detp)
optimize!(m_dccc_n2n_a_det)
termination_status(m_dccc_n2n_a_det)
objective_value(m_dccc_n2n_a_det)
value(m_dccc_n2n_a_det[:det_c])
value(m_dccc_n2n_a_det[:d_lin])
value(m_dccc_n2n_a_det[:d_con])
value(m_dccc_n2n_a_det[:d_quad])
value(m_dccc_n2n_a_det[:unc_c])
value(m_dccc_n2n_a_det[:u_quad])
value(m_dccc_n2n_a_det[:u_bil])

min_alg_a(10, α_min_initm, α_min_initp, lim_init = 0.00001, lim_const = 0.00000002)


# results_approx = Vector{Float64}()
# my_aprxs = Vector{aprx}()

# for i in 1:16
#
#     global my_aprxs = Vector{aprx}()
#     for j in 1:n_generators
#         push!(my_aprxs, aprx(approx(generators, j, segments_pu = i)...))
#     end
#
#     # To approximate the quadratic cost function, linear segments are used.
#     #     These are added as constraints. As these are greater than the function
#
#     include("models/dccc_sym_lin.jl")
#     m_dccc_lin = build_dccc_sym_lin(generators, buses, lines, farms)
#     optimize!(m_dccc_lin)
#     termination_status(m_dccc_lin)
#     push!(results_approx, objective_value(m_dccc_lin))
#
# end

#include("plot_approx.jl")

#-----------------------------------------------------
#-----------------------------------------------------
mm = Model()
@variable(mm, kp[1:n_generators] >= 0)
@variable(mm, km[1:n_generators] >= 0)

vec(vcat(kp, km)) * [μm..., μp]
[μm..., -μm...]


case_data, generators = updateGen(0.1, 0.39)

include("models/dccc.jl")
m_dccc = build_dccc(generators, buses, lines, farms)
optimize!(m_dccc)
termination_status(m_dccc)
z1 = objective_value(m_dccc)

include("models/dccc_sym.jl")
m_dccc = build_dccc_sym(generators, buses, lines, farms)
optimize!(m_dccc)
termination_status(m_dccc)
z2 = objective_value(m_dccc)

print(z1-z2)

for i in 1:11
    print(string(i, " & "))
end

for i in 1:11
    print(string(round(σ_vec[i], digits = 3), " & "))
end

for i in 1:11
    print(string(round(μm[i], digits = 3), " & "))
end

for i in 1:11
    print(string(round(σm[i], digits = 3), " & "))
end

print(sqrt(ζ2))

## System-Wide VS Node-To-Node
##############################

case_data, generators = updateGen(0.1, 0.36)

include("models/dccc_sym.jl")
m_dccc = build_dccc_sym(generators, buses, lines, farms)
optimize!(m_dccc)
termination_status(m_dccc)
z1 = objective_value(m_dccc)
a_s = value.(m_dccc[:α]) #* sum(σ_vec)
λ_s = -dual.(m_dccc[:mc])
γs = dual.(m_dccc[:γ])
p_s = value.(m_dccc[:p])
cc1 = -dual.(m_dccc[:cc1])
cc2 = -dual.(m_dccc[:cc2])
y = [z * value.(m_dccc[:α])[i] * s for i in 1:n_generators]
mi = y .+ [g.Pgmin for g in generators]
ma = [g.Pgmax for g in generators] .- y

include("models/dccc_n2n_sym.jl")
m_dccc_n2n = build_dccc_n2n_sym(generators, buses, lines, farms)
optimize!(m_dccc_n2n)
termination_status(m_dccc_n2n)
z2 = objective_value(m_dccc_n2n)
λ_s_n2n = -dual.(m_dccc_n2n[:mc])
α = value.(m_dccc_n2n[:α])
sum(value.(m_dccc_n2n[:p_uncert]))
χs = dual.(m_dccc_n2n[:χ])
sum(dual.(m_dccc_n2n[:χ]))

include("plot_sys_n2n.jl")

for i in χ
    print(string(round(i, digits = 2),"&"))
end
for i in χ
    print(string("-","&"))
end

## Symmetric VS Asymmetric
##########################

case_data, generators = updateGen(0.1, 0.39)

include("models/dccc.jl")
m_dccc = build_dccc(generators, buses, lines, farms)
optimize!(m_dccc)
termination_status(m_dccc)
z1 = objective_value(m_dccc)
a_s = value.(m_dccc[:α])
λ_s = -dual.(m_dccc[:mc])
γs = dual.(m_dccc[:γ])
p_s = value.(m_dccc[:p])
cc1 = -dual.(m_dccc[:cc1])
cc2 = -dual.(m_dccc[:cc2])
y = [z * value.(m_dccc[:α])[i] * s for i in 1:n_generators]

include("models/dccc_ab.jl")
m_dccc_ab = build_dccc_ab(generators, buses, lines, farms)
optimize!(m_dccc_ab)
termination_status(m_dccc_ab)
z3 = objective_value(m_dccc_ab)
sum(value.(m_dccc_ab[:ucp]))
sum(value.(m_dccc_ab[:cp]))
z3_up = value.(m_dccc_ab[:ucp])
z3_um = value.(m_dccc_ab[:cp])
ap = value.(m_dccc_ab[:αp])
am = value.(m_dccc_ab[:αm])
λ_ab  = -dual.(m_dccc_ab[:mc])
γp = dual.(m_dccc_ab[:γp])
γm = dual.(m_dccc_ab[:γm])
cc1 = dual.(m_dccc_ab[:cc1])
cc2 = dual.(m_dccc_ab[:cc2])

include("models/dccc_n2n.jl")
m_dccc_n2n = build_dccc_n2n(generators, buses, lines, farms)
optimize!(m_dccc_n2n)
termination_status(m_dccc_n2n)
z2 = objective_value(m_dccc_n2n)
λ_s_n2n = -dual.(m_dccc_n2n[:mc])
α = value.(m_dccc_n2n[:α])
p_u = value.(m_dccc_n2n[:p_uncert])
χs = dual.(m_dccc_n2n[:χ])
sum(dual.(m_dccc_n2n[:χ]))

include("models/dccc_n2n_ab.jl")
m_dccc_n2n_ab = build_dccc_n2n_ab(generators, buses, lines, farms)
optimize!(m_dccc_n2n_ab)
getobjectivevalue(m_dccc_n2n_ab)
termination_status(m_dccc_n2n_ab)
sum(value.(m_dccc_n2n_ab[:cp]))
sum(value.(m_dccc_n2n_ab[:ecp]))
z4 = objective_value(m_dccc_n2n_ab)
χp = dual.(m_dccc_n2n_ab[:χp])
sum(χp)
χm = dual.(m_dccc_n2n_ab[:χm])
sum(χm)
λ_n2n_ab = -dual.(m_dccc_n2n_ab[:mc])
pp_u = value.(m_dccc_n2n_ab[:pp_uncert])
pm_u = value.(m_dccc_n2n_ab[:pm_uncert])
value.(m_dccc_n2n_ab[:norm_up])
value.(m_dccc_n2n_ab[:norm_dwn])
value.(m_dccc_n2n_ab[:cc1])
value.(m_dccc_n2n_ab[:cc2])
am_n2n_ab = value.(m_dccc_n2n_ab[:αm])
ap_n2n_ab = value.(m_dccc_n2n_ab[:αp])
cc1_ab = dual.(m_dccc_ab[:cc1])
cc2_ab = dual.(m_dccc_ab[:cc2])

for i in 1:11
    print(string(round(χm[i], digits = 1)," & "))
end

for i in 1:11
    print(string(round(abs(χp[i]), digits = 1)," & "))
end

for i in 1:11
    print(string(0.0," & "))
end
# sum(sum(μm) * am + za * sm * am)
# sum(a_s * z * s) * ζ3

include("plot_sym_asym.jl")

## Symmetric VS Asymmetric
##########################

case_data, generators = updateGen(0.1, 0.365)

include("models/dccc.jl")
m_dccc = build_dccc(generators, buses, lines, farms)
optimize!(m_dccc)
termination_status(m_dccc)
z1 = objective_value(m_dccc)
a_s = value.(m_dccc[:α])
λ_s_c = -dual.(m_dccc[:mc])
γs = dual.(m_dccc[:γ])
p_s = value.(m_dccc[:p])
cc1s = -dual.(m_dccc[:cc1])
cc2s = -dual.(m_dccc[:cc2])
y = [z * value.(m_dccc[:α])[i] * s for i in 1:n_generators]

include("models/dccc_ab.jl")
m_dccc_ab = build_dccc_ab(generators, buses, lines, farms)
optimize!(m_dccc_ab)
termination_status(m_dccc_ab)
z3 = objective_value(m_dccc_ab)
sum(value.(m_dccc_ab[:ucp]))
sum(value.(m_dccc_ab[:cp]))
z3_up = value.(m_dccc_ab[:ucp])
z3_um = value.(m_dccc_ab[:cp])
ap = value.(m_dccc_ab[:αp])
am = value.(m_dccc_ab[:αm])
λ_ab  = -dual.(m_dccc_ab[:mc])
γp = dual.(m_dccc_ab[:γp])
γm = dual.(m_dccc_ab[:γm])
cc1s_ab = dual.(m_dccc_ab[:cc1])
cc2s_ab = dual.(m_dccc_ab[:cc2])

include("models/dccc_n2n.jl")
m_dccc_n2n = build_dccc_n2n(generators, buses, lines, farms)
optimize!(m_dccc_n2n)
termination_status(m_dccc_n2n)
z2 = objective_value(m_dccc_n2n)
λ_s_n2n = -dual.(m_dccc_n2n[:mc])
α = value.(m_dccc_n2n[:α])
p_u = value.(m_dccc_n2n[:p_uncert])
χs = dual.(m_dccc_n2n[:χ])
sum(dual.(m_dccc_n2n[:χ]))
cc1 = dual.(m_dccc_n2n[:cc1])
cc2 = dual.(m_dccc_n2n[:cc2])

include("models/dccc_n2n_ab.jl")
m_dccc_n2n_ab = build_dccc_n2n_ab(generators, buses, lines, farms)
optimize!(m_dccc_n2n_ab)
termination_status(m_dccc_n2n_ab)
sum(value.(m_dccc_n2n_ab[:cp]))
sum(value.(m_dccc_n2n_ab[:ecp]))
z4 = objective_value(m_dccc_n2n_ab)
χp = dual.(m_dccc_n2n_ab[:χp])
sum(χp)
χm = dual.(m_dccc_n2n_ab[:χm])
sum(χm)
λ_n2n_ab_c = -dual.(m_dccc_n2n_ab[:mc])
pp_u = value.(m_dccc_n2n_ab[:pp_uncert])
pm_u = value.(m_dccc_n2n_ab[:pm_uncert])
value.(m_dccc_n2n_ab[:norm_up])
value.(m_dccc_n2n_ab[:norm_dwn])
value.(m_dccc_n2n_ab[:cc1])
value.(m_dccc_n2n_ab[:cc2])
am_n2n_ab = value.(m_dccc_n2n_ab[:αm])
ap_n2n_ab = value.(m_dccc_n2n_ab[:αp])
cc1_ab_c = dual.(m_dccc_n2n_ab[:cc1])
cc2_ab_c = dual.(m_dccc_n2n_ab[:cc2])

include("save_data.jl")
include("plot_delta.jl")
