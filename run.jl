cd(@__DIR__)
include("pkgs.jl")
include("code_jl/input.jl")
include("code_jl/gradients.jl")
#include("code_jl/TruncatedGaussian.jl")

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

sum([g.Pgmax for g in generators])

# for i in σ_vec
#     print(string(i,", "))
# end

B = X ^ (-1) * A
B_node = A' * B

d = [b.Pd for b in buses]
sum([b.Pd for b in buses])
sum([g.Pgmax for g in generators])

# include("code_jl/linApprox.jl")
using JLD
using Distributions
include("code_jl//farms.jl")
uRESs, n_ures, σ, Σ, s_sq, Σ_rt, s, μ = create_wind_farms(fc = 2)

sum([u.forecast for u in uRESs])

μs = sum(μ)
𝛭 = sum(μ[f] for f in 1:length(uRESs))

ζ1 = (sqrt(pi / 2))
ζ2 = (2 * pi - 4) / (2 * pi)

truncated_dists_p = Vector{Distribution}()
truncated_dists_m = Vector{Distribution}()
for i in 1:length(uRESs)
    push!(truncated_dists_p, TruncatedNormal(μ[i], σ[i], 0, Inf64))
    push!(truncated_dists_m, TruncatedNormal(μ[i], σ[i], -Inf64, 0))
end

u = 4
mean(truncated_dists_p[u])
mean(truncated_dists_m[u])
sqrt(var(truncated_dists_p[u]))
sqrt(var(truncated_dists_m[u]))

σm = [sqrt(var(d)) for d in truncated_dists_m]
σp = [sqrt(var(d)) for d in truncated_dists_p]

μm = [mean(d) for d in truncated_dists_m] .* (-1)
μp = [mean(d) for d in truncated_dists_p]

# σm = σ .* sqrt(ζ2)
# σp = σ .* sqrt(ζ2)
#
# μm = μ .+ (σ .* (1 / ζ1))
# μp = μ .- (σ .* (1 / ζ1))

sm = sqrt(sum(σm.^2))
sp = sqrt(sum(σp.^2))

μms = sum(μm)
μps = sum(μp)

Σ = (diagm(σ))^2
Σm = diagm(σm).^2
Σp = diagm(σp).^2

Σ_rt = sqrt(Σ)
Σm_rt = sqrt(Σm)
Σp_rt = sqrt(Σp)

n_farms = length(uRESs)

#include("plotDists.jl")
# v = 7
# println(z * σ[v])
# println(μm[v] + za * σm[v])
# println("----------")

d = [b.Pd for b in buses]

counter = 1
for u in uRESs
    print(string(round(u.forecast, digits = 3)," & "))
end

println("\n mu \n")

counter = 1
for u in uRESs
    print(string(round(u.μ, digits = 3)," & "))
end

println("\n mum \n")

counter = 1
for u in uRESs
    print(string(round(μm[counter], digits = 3)," & "))
    global counter += 1
end

println("\n mup \n")

counter = 1
for u in uRESs
    print(string(round(μp[counter], digits = 3)," & "))
    global counter += 1
end

println("\n sigma \n")

counter = 1
for u in uRESs
    print(string(round(u.σ, digits = 3)," & "))
end

println("\n sigmam \n")

counter = 1
for u in uRESs
    print(string(round(σm[counter], digits = 3)," & "))
    global counter += 1
end

println("\n sigmap \n")

counter = 1
for u in uRESs
    print(string(round(σp[counter], digits = 3)," & "))
    global counter += 1
end



#
# counter = 1
# for f in farms
#     print(string(round(f.σ, digits = 4)," & "))
#     global counter += 1
# end
#
# counter = 1
# for f in farms
#     print(string(round(μm[counter], digits = 4)," & "))
#     global counter += 1
# end
#
# counter = 1
# for f in farms
#     print(string(round(Σm_rt[counter,counter], digits = 4)," & "))
#     global counter += 1
# end
#
# for (i,f) in enumerate(farms)
#     push!(buses[f.bus].farmids, i)
# end

## Stochastic parameters
########################

ϵ = 0.01
z = quantile(Normal(0,1), 1-ϵ)
za = quantile(Normal(0,1), 1 - ϵ / 0.5)
σ_cheb = sqrt((1 - ϵ)/ϵ) .* σ .- μ

## Generation Costs
###################

c_vec = [g.pi1  for g in generators]
C_mat = diagm(0 => c_vec)
C_rt = sqrt(C_mat)

#include("Approx_planes.jl")

@inline function updateGen(min::Float64, max::Float64)

    @assert min >= 0 && min <= 1 && max <= 1 && max >= 0 "Values must be between 0 and 1."

    case_data = load("data//118bus.jld")
    generators = case_data["generators"]

    for (i, g) in enumerate(generators)
         g.Pgmin = g.Pgmax * min
         g.Pgmax = g.Pgmax * max
    end

    global generators = generators
    #return case_data, generators
end

@inline function scaleCosts(factor::Float64)

    for (i, g) in enumerate(generators)
         g.pi1 = g.pi1 / sqrt(factor)
         g.pi2 = g.pi2 / factor
         g.pi3 = g.pi3 / factor
    end

    global generators = generators
    #return case_data, generators
end

## Experiments
##############

## Linear Costs
###############
updateGen(0.1, 0.6)
scaleCosts(1.0)

include("models/dccc.jl")
m_dccc = build_dccc(generators, buses, lines, uRESs)
optimize!(m_dccc)
termination_status(m_dccc)
z1 = objective_value(m_dccc)
value(m_dccc[:det_c])
value(m_dccc[:d_lin])
value(m_dccc[:d_con])
value(m_dccc[:d_quad])
value(m_dccc[:unc_c])
χ0 = dual.(m_dccc[:γ])
δ1ms = dual.(m_dccc[:cc1])
δ1ps = dual.(m_dccc[:cc2])

include("models/dccc_apx.jl")
mylim = 0.0000
m_dccc_apx = build_dccc_apx(generators, buses, lines, uRESs, α_min = [mylim for i in 1:n_generators])
optimize!(m_dccc_apx)
termination_status(m_dccc_apx)
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
χ1 = dual.(m_dccc_apx[:γ])

include("models/dccc_det.jl")
m_dccc_det = build_dccc_det(generators, buses, lines, uRESs, value.(m_dccc_apx[:p]))
optimize!(m_dccc_det)
termination_status(m_dccc_det)
z11 = objective_value(m_dccc_det)
value(m_dccc_det[:det_c])
value(m_dccc_det[:d_lin])
value(m_dccc_det[:d_con])
value(m_dccc_det[:d_quad])
value(m_dccc_det[:unc_c])
value(m_dccc_det[:u_bil])
value(m_dccc_det[:u_quads])
value(m_dccc_det[:u_quadm])
χ1 = dual.(m_dccc_det[:γ])
δ1ms = dual.(m_dccc_det[:cc1])
δ1ps = dual.(m_dccc_det[:cc2])
sum(δ1ms)
sum(δ1ps)

include("models/dccc_a_apx.jl")
m_dccc_a_apx = build_dccc_a_apx(generators, buses, lines, uRESs, αm_min = zeros(n_generators), αp_min = zeros(n_generators))
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
sum(value.(m_dccc_a_apx[:p]))
sum(value.(m_dccc_det[:p]))
value.(m_dccc_a_apx[:p])
d2 = dual.(m_dccc_a_apx[:γm])
d2 = dual.(m_dccc_a_apx[:γp])
value.(m_dccc_a_apx[:p])

include("models/dccc_a_det.jl")
m_dccc_a_det = build_dccc_a_det(generators, buses, lines, uRESs, value.(m_dccc_a_apx[:p]))
optimize!(m_dccc_a_det)
termination_status(m_dccc_a_det)
z4 = objective_value(m_dccc_a_det)
value(m_dccc_a_det[:det_c])
value(m_dccc_a_det[:d_lin])
value(m_dccc_a_det[:d_con])
value(m_dccc_a_det[:d_quad])
value(m_dccc_a_det[:u_quads_m])
value(m_dccc_a_det[:u_quads_p])
value(m_dccc_a_det[:u_quadm_m])
value(m_dccc_a_det[:u_quadm_p])
value(m_dccc_a_det[:unc_c])
value(m_dccc_a_det[:d_bil])
χm = dual.(m_dccc_a_det[:γm])
χp = dual.(m_dccc_a_det[:γp])

value.(m_dccc_det[:α])
value.(m_dccc_a_det[:αp])

function sw2n2n(alpha)

    # g x f
    new_alpha = ones(n_generators, length(uRESs))
    for i in 1:length(uRESs)
        new_alpha[:, i] = [alpha[i] for i in 1:n_generators]
    end
    return new_alpha

end

include("models/dccc_n2n.jl")
m_dccc_n2n = build_dccc_n2n(generators, buses, lines, uRESs)
optimize!(m_dccc_n2n)
termination_status(m_dccc_n2n)
z4 = objective_value(m_dccc_n2n)
value(m_dccc_n2n[:det_c])
value(m_dccc_n2n[:d_lin])
value(m_dccc_n2n[:d_con])
value(m_dccc_n2n[:d_quad])
value(m_dccc_n2n[:unc_c])
sum(value.(m_dccc_n2n[:p_uncert]) .^ 2)
χ0u = dual.(m_dccc_n2n[:χ])
δ0m = dual.(m_dccc_n2n[:cc1])
δ0p = dual.(m_dccc_n2n[:cc2])
sum(δ0m)
sum(δ0p)


sum(χ1u)
for i in 1:11
    print(string(round(χ1u[i], digits = 1)," & "))
end

α_min_init = zeros((n_generators, length(uRESs))) .+ 0.003 # sw2n2n(value.(m_dccc_a_det[:αp])) .- 0.000 #ones((n_generators, n_farms)) * 0.000001
α_max_init = ones((n_generators, length(uRESs)))
include("models/dccc_n2n_apx.jl")
m_dccc_n2n_apx = build_dccc_n2n_apx(generators, buses, lines, uRESs, α_min_init, α_max_init)
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
d31 = dual.(m_dccc_n2n_apx[:χ])

α_det = value.(m_dccc_n2n_apx[:α])
p_det = value.(m_dccc_n2n[:p])

include("models/dccc_n2n_det.jl")
m_dccc_n2n_det = build_dccc_n2n_det(generators, buses, lines, uRESs, α_det, p_det = p_det)
optimize!(m_dccc_n2n_det )
termination_status(m_dccc_n2n_det)
z4 = objective_value(m_dccc_n2n_det)
value(m_dccc_n2n_det[:det_c])
value(m_dccc_n2n_det[:d_lin])
value(m_dccc_n2n_det[:d_con])
value(m_dccc_n2n_det[:d_quad])
value(m_dccc_n2n_det[:unc_c])
value(m_dccc_n2n_det[:u_bil])
value(m_dccc_n2n_det[:u_quads])
value(m_dccc_n2n_det[:u_quadm])
χ1u = dual.(m_dccc_n2n_det[:χ])
δ1m = dual.(m_dccc_n2n_det[:cc1])
δ1p = dual.(m_dccc_n2n_det[:cc2])
sum(δ1m)
sum(δ1p)

α_min_initm = zeros((n_generators, n_farms)) .+ 0.0 # value.(m_dccc_n2n_det[:α]) # ones((n_generators, n_farms)) * 0.01
α_max_initm = ones((n_generators, n_farms)) * 1.0
α_min_initp = zeros((n_generators, n_farms)) .+ 0.0# value.(m_dccc_n2n_det[:α]) # ones((n_generators, n_farms)) * 0.01
α_max_initp = ones((n_generators, n_farms)) * 1.0

include("models/dccc_n2n_a_apx_alpha.jl")
m_dccc_n2n_a_apx_alpha = build_dccc_n2n_a_apx_alpha(generators, buses, lines, uRESs, α_min_initm, α_max_initm, α_min_initp, α_max_initp)
optimize!(m_dccc_n2n_a_apx_alpha)
objective_value(m_dccc_n2n_a_apx_alpha)
termination_status(m_dccc_n2n_a_apx_alpha)
value(m_dccc_n2n_a_apx_alpha[:det_c])
value(m_dccc_n2n_a_apx_alpha[:d_lin])
value(m_dccc_n2n_a_apx_alpha[:d_con])
value(m_dccc_n2n_a_apx_alpha[:d_quad])
value(m_dccc_n2n_a_apx_alpha[:unc_c])
value(m_dccc_n2n_a_apx_alpha[:u_bil])
d41 = dual.(m_dccc_n2n_a_apx_alpha[:χp])
d42 = dual.(m_dccc_n2n_a_apx_alpha[:χm])

α_detm = value.(m_dccc_n2n_a_apx_alpha[:αm])
α_detp = value.(m_dccc_n2n_a_apx_alpha[:αp])
p_det = value.(m_dccc_n2n_a_apx_alpha[:p])
# p_det = value.(m_dccc_n2n_det[:p])

include("models/dccc_n2n_a_det_p.jl")
m_dccc_n2n_a_det_p = build_dccc_n2n_a_det_p(generators, buses, lines, uRESs, p_det)
optimize!(m_dccc_n2n_a_det_p)
termination_status(m_dccc_n2n_a_det_p)
z4 = objective_value(m_dccc_n2n_a_det_p)
value(m_dccc_n2n_a_det_p[:det_c])
value(m_dccc_n2n_a_det_p[:d_lin])
value(m_dccc_n2n_a_det_p[:d_con])
value(m_dccc_n2n_a_det_p[:d_quad])
value(m_dccc_n2n_a_det_p[:unc_c])
value(m_dccc_n2n_a_det_p[:u_quad])
value(m_dccc_n2n_a_det_p[:u_bil])
χ1up = dual.(m_dccc_n2n_a_det_p[:χp])
χ1um = dual.(m_dccc_n2n_a_det_p[:χm])

include("models/dccc_n2n_a_det_alpha.jl")

# min_alg_a(10, α_min_initm, α_min_initp, lim_init = 0.00001, lim_const = 0.00000002)
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

include("plot_sym_asym1.jl")

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
