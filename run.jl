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
σ
ϵ = 0.01
z = quantile(Normal(0,1), 1-ϵ)
za = quantile(Normal(0,1), 1 - ϵ / 0.5)
σ_cheb = sqrt((1 - ϵ)/ϵ) .* σ .- μ
σ_cheb_m = sqrt((1 - ϵ)/ϵ) .* σm .- μm
σ_cheb_p = sqrt((1 - ϵ)/ϵ) .* σp .- μp

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

# Parameters
# Dual values
σ_cheb
σ_cheb_m
σ_cheb_p

counter = 1
for f in uRESs
    print(string(round(σ_cheb[counter], digits = 4)," & "))
    global counter += 1
end

counter = 1
for f in uRESs
    print(string(round(σ_cheb_m[counter], digits = 4)," & "))
    global counter += 1
end

counter = 1
for f in uRESs
    print(string(round(σ_cheb_p[counter], digits = 4)," & "))
    global counter += 1
end

function get_payments(m)

    generation = zeros(length(buses))
    genids = [g.busidx for g in generators]

    ctr = 1
    for g in genids
        generation[g] = value.(m[:p])[ctr]
        ctr += 1
    end

    generation_fc = zeros(length(buses))

    sum(dual.(m[:mc]) .* d * (-1)) - sum(generation_fc .* dual.(m[:mc]) * (-1))
    sum(-generation .* dual.(m[:mc]))

    ctr = 1
    for g in [ur.bus for ur in uRESs]
        generation_fc[g] = [ur.forecast for ur in uRESs][ctr]
        ctr += 1
    end

           # consumer payments,             wind payments,                              generator payments
    return sum(dual.(m[:mc]) .* d * (-1)), sum(generation_fc .* dual.(m[:mc]) * (-1)), sum(-generation .* dual.(m[:mc]))
end

include("models/dc.jl")
m_dc = build_dc(generators, buses, lines, uRESs)
optimize!(m_dc)
termination_status(m_dc)
z_dc = objective_value(m_dc)
value(m_dc[:det_c])
value(m_dc[:d_lin])
value(m_dc[:d_con])
value(m_dc[:d_quad])
χ0 = dual.(m_dc[:γ])
δ1ms = dual.(m_dc[:cc1])
δ1ps = dual.(m_dc[:cc2])

# pdet = value.(m_dc[:p])
# adet = value.(m_dc[:α])
#
# zdcp = get_payments(m_dc)
# zdcp2 = get_payments(m_dc)
# zdcp3 = get_payments(m_dc)

include("models/dccc.jl")
m_dccc = build_dccc(generators, buses, lines, uRESs)
optimize!(m_dccc)
termination_status(m_dccc)
z1n = objective_value(m_dccc)
z1ns = objective_value(m_dccc) - value(m_dccc[:unc_c])
value(m_dccc[:det_c])
value(m_dccc[:d_lin])
value(m_dccc[:d_con])
value(m_dccc[:d_quad])
value(m_dccc[:unc_c])
χ0 = dual.(m_dccc[:γ])
δ1ms = dual.(m_dccc[:cc1])
δ1ps = dual.(m_dccc[:cc2])

z1np = get_payments(m_dccc_det)

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
z1 = objective_value(m_dccc_det)
z1s = objective_value(m_dccc_det) - value(m_dccc_det[:unc_c])
value(m_dccc_det[:det_c])
value(m_dccc_det[:d_lin])
value(m_dccc_det[:d_con])
value(m_dccc_det[:d_quad])
value(m_dccc_det[:unc_c])
value(m_dccc_det[:u_quads])
value(m_dccc_det[:u_quadm])
χ1 = dual.(m_dccc_det[:γ])
δ1ms = dual.(m_dccc_det[:cc1])
δ1ps = dual.(m_dccc_det[:cc2])
sum(δ1ms)
sum(δ1ps)

z1p = get_payments(m_dccc_det)

#updateGen(0.1, 0.6)
include("models/dccc_det_cheb.jl")
m_dccc_det = build_dccc_det_cheb(generators, buses, lines, uRESs, value.(m_dccc_apx[:p]))
optimize!(m_dccc_det)
termination_status(m_dccc_det)
z1c = objective_value(m_dccc_det)
z1c_d = objective_value(m_dccc_det) - value(m_dccc_det[:unc_c])
value(m_dccc_det[:det_c])
value(m_dccc_det[:d_lin])
value(m_dccc_det[:d_con])
value(m_dccc_det[:d_quad])
value(m_dccc_det[:unc_c])
#value(m_dccc_det[:u_bil])
# value(m_dccc_det[:u_quads])
#value(m_dccc_det[:u_quadm])
χ1 = dual.(m_dccc_det[:γ])
δ1ms = dual.(m_dccc_det[:cc1])
δ1ps = dual.(m_dccc_det[:cc2])
sum(δ1ms)
sum(δ1ps)

z1cp = get_payments(m_dccc_det)

include("models/dccc_a_apx.jl")
m_dccc_a_apx = build_dccc_a_apx(generators, buses, lines, uRESs, αm_min = zeros(n_generators), αp_min = zeros(n_generators))
optimize!(m_dccc_a_apx)
termination_status(m_dccc_a_apx)
z2apx = objective_value(m_dccc_a_apx)
value(m_dccc_a_apx[:det_c])
value(m_dccc_a_apx[:d_lin])
value(m_dccc_a_apx[:d_con])
value(m_dccc_a_apx[:d_quad])
value(m_dccc_a_apx[:unc_c])
value(m_dccc_a_apx[:d_bil])
value.(m_dccc_a_apx[:αm])
value.(m_dccc_a_apx[:αp])
sum(value.(m_dccc_a_apx[:p]))
# sum(value.(m_dccc_det[:p]))
value.(m_dccc_a_apx[:p])
d2 = dual.(m_dccc_a_apx[:γm])
d2 = dual.(m_dccc_a_apx[:γp])
value.(m_dccc_a_apx[:p])

apdet = value.(m_dccc_a_apx[:αm])
amdet = value.(m_dccc_a_apx[:αp])

include("models/dccc_a_det.jl")
#m_dccc_a_det = build_dccc_a_det(generators, buses, lines, uRESs; pdet = value.(m_dccc_a_apx[:p]))
m_dccc_a_det = build_dccc_a_det(generators, buses, lines, uRESs; apdet = apdet, amdet = amdet)
# m_dccc_a_det = build_dccc_a_det(generators, buses, lines, uRESs)
optimize!(m_dccc_a_det)
termination_status(m_dccc_a_det)
z2 = objective_value(m_dccc_a_det)
z2s = objective_value(m_dccc_a_det) - value(m_dccc_a_det[:unc_c])
value(m_dccc_a_det[:det_c])
value(m_dccc_a_det[:d_lin])
value(m_dccc_a_det[:d_con])
value(m_dccc_a_det[:d_quad])
value(m_dccc_a_det[:u_quads_m])
value(m_dccc_a_det[:u_quads_p])
value(m_dccc_a_det[:u_quadm_m])
value(m_dccc_a_det[:u_quadm_p])
value(m_dccc_a_det[:unc_c])
# value(m_dccc_a_det[:d_bil])
# χm = dual.(m_dccc_a_det[:γm])
# χp = dual.(m_dccc_a_det[:γp])

dual_objective_value(m_dccc_a_det)

m_dccc_a_det[:eq22]

sum(dual.(m_dccc_a_det[:mc]))
apdet = abs.(value.(m_dccc_a_det[:αp]))
amdet = abs.(value.(m_dccc_a_det[:αm]))

include("models/dccc_a_det_falpha.jl")
m_dccc_a_det_fa = build_dccc_a_det_falpha(generators, buses, lines, uRESs; apdet = apdet, amdet = amdet)
optimize!(m_dccc_a_det_fa)
objective_value(m_dccc_a_det_fa)
termination_status(m_dccc_a_det_fa)
dual.(m_dccc_a_det_fa[:mc])
value.(m_dccc_a_det_fa[:p])

z2p = get_payments(m_dccc_a_det_fa)


dual.(m_dccc_a_det[:mc])
z2p = get_payments(m_dccc_a_det)
z2p2 = get_payments(m_dccc_a_det)

generation = zeros(length(buses))
genids = [g.busidx for g in generators]

ctr = 1
for g in genids
    generation[g] = value.(m_dccc_a_det[:p])[ctr]
    global ctr += 1
end

dual.(m_dccc_a_det[:mc]) .* generation

function get_payments(m)

    generation = zeros(length(buses))
    genids = [g.busidx for g in generators]

    ctr = 1
    for g in genids
        generation[g] = value.(m[:p])[ctr]
        ctr += 1
    end

    generation_fc = zeros(length(buses))

    sum(dual.(m[:mc]) .* d * (-1)) - sum(generation_fc .* dual.(m[:mc]) * (-1))
    sum(-generation .* dual.(m[:mc]))

    ctr = 1
    for g in [ur.bus for ur in uRESs]
        generation_fc[g] = [ur.forecast for ur in uRESs][ctr]
        ctr += 1
    end

           # consumer payments,             wind payments,                              generator payments
    return sum(dual.(m[:mc]) .* d * (-1)), sum(generation_fc .* dual.(m[:mc]) * (-1)), sum(-generation .* dual.(m[:mc]))
end

include("models/dccc_a_det_cheb.jl")
m_dccc_a_det_cheb = build_dccc_a_det_cheb(generators, buses, lines, uRESs, value.(m_dccc_a_apx[:p]))
optimize!(m_dccc_a_det_cheb)
termination_status(m_dccc_a_det_cheb)
z2c = objective_value(m_dccc_a_det_cheb)
z2c_d = objective_value(m_dccc_a_det_cheb) - value(m_dccc_a_det_cheb[:unc_c])
value(m_dccc_a_det_cheb[:det_c])
value(m_dccc_a_det_cheb[:d_lin])
value(m_dccc_a_det_cheb[:d_con])
value(m_dccc_a_det_cheb[:d_quad])
value(m_dccc_a_det_cheb[:u_quads_m])
value(m_dccc_a_det_cheb[:u_quads_p])
# value(m_dccc_a_det_cheb[:u_quadm_m])
# value(m_dccc_a_det_cheb[:u_quadm_p])
value(m_dccc_a_det_cheb[:unc_c])
#value(m_dccc_a_det_cheb[:d_bil])
χm = dual.(m_dccc_a_det_cheb[:γm])
χp = dual.(m_dccc_a_det_cheb[:γp])

value.(m_dccc_det[:α])
value.(m_dccc_a_det[:αp])

z2cp = get_payments(m_dccc_a_det_cheb)

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
z3n = objective_value(m_dccc_n2n)
z3ns = objective_value(m_dccc_n2n) - value(m_dccc_n2n[:unc_c])
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

z3np = get_payments(m_dccc_n2n)
z1np
# sum(χ1u)
# for i in 1:11
#     print(string(round(χ1u[i], digits = 1)," & "))
# end

α_min_init = zeros((n_generators, length(uRESs))) .+ 0.000 # sw2n2n(value.(m_dccc_a_det[:αp])) .- 0.000 #ones((n_generators, n_farms)) * 0.000001
α_max_init = ones((n_generators, length(uRESs)))
include("models/dccc_n2n_apx.jl")
m_dccc_n2n_apx = build_dccc_n2n_apx(generators, buses, lines, uRESs, α_min_init, α_max_init)
optimize!(m_dccc_n2n_apx)
termination_status(m_dccc_n2n_apx)
z3ap = objective_value(m_dccc_n2n_apx)
value(m_dccc_n2n_apx[:det_c])
value(m_dccc_n2n_apx[:d_lin])
value(m_dccc_n2n_apx[:d_con])
value(m_dccc_n2n_apx[:d_quad])
value(m_dccc_n2n_apx[:unc_c])
value(m_dccc_n2n_apx[:u_bil])
value.(m_dccc_n2n_apx[:p_uncert])
d31 = dual.(m_dccc_n2n_apx[:χ])

α_det = value.(m_dccc_n2n_apx[:α])
p_det = value.(m_dccc_n2n_apx[:p])

include("models/dccc_n2n_det.jl")
m_dccc_n2n_det = build_dccc_n2n_det(generators, buses, lines, uRESs, α_det, p_det = p_det)
optimize!(m_dccc_n2n_det)
termination_status(m_dccc_n2n_det)
z3 = objective_value(m_dccc_n2n_det)
z3s = objective_value(m_dccc_n2n_det) - value(m_dccc_n2n_det[:unc_c])
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

z3p = get_payments(m_dccc_n2n_det)

include("models/dccc_n2n_det_cheb.jl")
m_dccc_n2n_det_cheb = build_dccc_n2n_det_cheb(generators, buses, lines, uRESs, α_det, p_det = p_det)
optimize!(m_dccc_n2n_det_cheb)
termination_status(m_dccc_n2n_det_cheb)
z3c = objective_value(m_dccc_n2n_det_cheb)
z3c_d = objective_value(m_dccc_n2n_det_cheb) - value(m_dccc_n2n_det_cheb[:unc_c])
value(m_dccc_n2n_det_cheb[:det_c])
value(m_dccc_n2n_det_cheb[:d_lin])
value(m_dccc_n2n_det_cheb[:d_con])
value(m_dccc_n2n_det_cheb[:d_quad])
value.(m_dccc_n2n_det_cheb[:p_uncert])
value(m_dccc_n2n_det_cheb[:unc_c])
χ1uc = dual.(m_dccc_n2n_det_cheb[:χ])
δ1m = dual.(m_dccc_n2n_det_cheb[:cc1])
δ1p = dual.(m_dccc_n2n_det_cheb[:cc2])
sum(δ1m)
sum(δ1p)

z3cp = get_payments(m_dccc_n2n_det_cheb)

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
z4s = objective_value(m_dccc_n2n_a_det_p) - value(m_dccc_n2n_a_det_p[:unc_c])
value(m_dccc_n2n_a_det_p[:det_c])
value(m_dccc_n2n_a_det_p[:d_lin])
value(m_dccc_n2n_a_det_p[:d_con])
value(m_dccc_n2n_a_det_p[:d_quad])
value(m_dccc_n2n_a_det_p[:unc_c])
value(m_dccc_n2n_a_det_p[:u_quad])
value(m_dccc_n2n_a_det_p[:u_bil])
χ1up = dual.(m_dccc_n2n_a_det_p[:χp])
χ1um = dual.(m_dccc_n2n_a_det_p[:χm])

χ1up - χ1um

include("models/dccc_n2n_a_det_p_cheb.jl")
m_dccc_n2n_a_det_p_cheb = build_dccc_n2n_a_det_p_cheb(generators, buses, lines, uRESs)
optimize!(m_dccc_n2n_a_det_p_cheb)
termination_status(m_dccc_n2n_a_det_p_cheb)
z4c = objective_value(m_dccc_n2n_a_det_p_cheb)
z4c_d = objective_value(m_dccc_n2n_a_det_p_cheb) - value(m_dccc_n2n_a_det_p_cheb[:unc_c])
value(m_dccc_n2n_a_det_p_cheb[:det_c])
value(m_dccc_n2n_a_det_p_cheb[:d_lin])
value(m_dccc_n2n_a_det_p_cheb[:d_con])
value(m_dccc_n2n_a_det_p_cheb[:d_quad])
value(m_dccc_n2n_a_det_p_cheb[:unc_c])
χ1up_cheb = dual.(m_dccc_n2n_a_det_p_cheb[:χp])
χ1um_cheb = dual.(m_dccc_n2n_a_det_p_cheb[:χm])

z4cp = get_payments(m_dccc_n2n_a_det_p_cheb)

sched_gen = sum(value.(m_dccc_n2n_a_det_p_cheb[:p]))
ures_gen = sum([my_u.forecast for my_u in uRESs])
total_load = sum(d)
ures_pen = ures_gen / total_load

using PyPlot

fig = figure(figsize=(2.8, 2.0))
rc("font", family = "serif", style = "italic", size = 12)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.20,0.23,0.78,0.76])
grid(linewidth = 0.1, linestyle = (0, (20, 20)), color = "lightgray", axis = "y")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4,length = 2.0)

xlabel("Balancing Framework")
ylim(bottom = 40, top = 260)
#ylim(bottom = 65, top = 100)
#ylim(bottom = 65, top = 128)
#xlim(left = 1, right = n_generators)
ylabel("\$z^{*} [1000\\\$]\$")
yticks([60,80,100,120,140,160], labels = ["60","80","100","120","140","160"])
xticks([2,6,10,14], labels = ["S-SW","A-SW","S-N2N","A-N2N"])
#yticks([120,125,130], labels = ["120","125","130"])
#yticks([110,115,120], labels = ["110","115","120"])
#yticks([105,110,115,120,125], labels = ["105","110","115","120","125"])
#yticks([75,80,85], labels = ["75","80","85"])

#bar([3,7,11,15], [z1c_d, z2c_d, z3c_d, z4c_d] ./ 1000, color = "mediumblue", width = 0.6, label = "Chebychev, scheduled")
#bar([1,9], [z1ns, z3ns] ./ 1000, color = "silver", width = 0.6, label = "SOC, scheduled")
#bar([2,6,10,14], [z1s, z2s, z3s, z4s] ./ 1000, color = "coral", width = 0.6, label = "Non-Convex, scheduled")
#bar([1,2,3,6,7,9,10,11,14,15], [z1n, z1, z1c, z2, z2c, z3n, z3, z3c,z4,z4c] ./ 1000, color = "black", width = 0.15, label = "Obj. with balancing")

bar([1,5,9,13], [z1c, z2c, z3c, z4c] ./ 1000, color = [70,80,190] ./ 255, width = 1.0, label = "\$z^{*}\$")
bar([2,6,10,14], [z1cp[1], z2cp[1], z3cp[1], z4cp[1]] ./ 1000, color = [122,1,119] ./ 255, width = 1.0, label = "\$\\lambda D\$")
bar([3,7,11,15], [z1cp[2] + z1cp[3], z2cp[2] + z2cp[3], z3cp[2] + z3cp[3], z4cp[2] + z4cp[3]] ./ 1000, color = [15,170,160] ./ 255, width = 1.0,label = "\$\\lambda W\$" )
bar([3,7,11,15], [z1cp[3], z2cp[3], z3cp[3], z4cp[3]] ./ 1000, color = [40,190,90] ./ 255, width = 1.0,label = "\$\\sum_{i}\\Pi_{i}\$" )
annotate("infeas.", [2.0 - 0.6, 52], rotation = 90, color = "black")
annotate("infeas.", [6.0 - 0.6, 52], rotation = 90, color = "black")

legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 2)
savefig(string("plots_final//payments$(ures_pen).pdf"), format = :pdf)

1==1
# fig = figure(figsize=(5, 2.6))
# rc("font", family = "serif", style = "italic", size = 14)
# rc("text", usetex = true)
# rc("lines", linewidth = 1)
#
# ax = fig.add_axes([0.13,0.18,0.86,0.8])
# grid(linewidth = 0.1, linestyle = (0, (10, 10)), color = "lightgray")
# ax.tick_params(direction = "in", top = true, right = true, width = 1.4, length=0)
#
# xlabel("Balancing Framework")
# ylim(bottom = 95, top = 139.5)
# #xlim(left = 1, right = n_generators)
# ylabel("\$z^{*} [1000\\\$]\$")
# xticks([2,6,10,14], labels = ["S-SW","A-SW","S-N2N","A-N2N"])
#
# bar([3,7,11,15], [z1c_d, z2c_d, z3c_d, z4c_d] ./ 1000, color = "mediumblue", width = 0.6, label = "Chebychev, scheduled")
# bar([1,9], [z1ns, z3ns] ./ 1000, color = "silver", width = 0.6, label = "SOC, scheduled")
# bar([2,6,10,14], [z1s, z2s, z3s, z4s] ./ 1000, color = "coral", width = 0.6, label = "Non-Convex, scheduled")
# bar([1,2,3,6,7,9,10,11,14,15], [z1n, z1, z1c, z2, z2c, z3n, z3, z3c,z4,z4c] ./ 1000, color = "black", width = 0.15, label = "Obj. with balancing")
#
# legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9)
# savefig(string("plots_final//objectives4.pdf"), format = :pdf)
#
# ##----------------------
# ##----------------------
#
# # [240,249,232] ./ 255
# # 204,235,197
# # 168,221,181
# # 123,204,196
# # 67,162,202
# # 8,104,172
#
# fig = figure(figsize=(6, 2.6))
# rc("font", family = "serif", style = "italic", size = 14)
# rc("text", usetex = true)
# rc("lines", linewidth = 1)
#
# ax = fig.add_axes([0.10,0.18,0.89,0.8])
# grid(linewidth = 0.1, linestyle = (0, (20, 20)), color = "lightgray", axis = "y")
# ax.tick_params(direction = "in", top = true, right = true, width = 1.4,length = 2.0)
#
# xlabel("Balancing Framework")
# ylim(bottom = 120, top = 160)
# #ylim(bottom = 65, top = 128)
# #xlim(left = 1, right = n_generators)
# ylabel("\$z^{*} [1000\\\$]\$")
# xticks([2.0,6.5,12.5,18.5,24.5], labels = ["DC", "S-SW","A-SW","S-N2N","A-N2N"])
# yticks([130,135,140], labels = ["130","135","140"])
#
# #bar([3,7,11,15], [z1c_d, z2c_d, z3c_d, z4c_d] ./ 1000, color = "mediumblue", width = 0.6, label = "Chebychev, scheduled")
# #bar([1,9], [z1ns, z3ns] ./ 1000, color = "silver", width = 0.6, label = "SOC, scheduled")
# #bar([2,6,10,14], [z1s, z2s, z3s, z4s] ./ 1000, color = "coral", width = 0.6, label = "Non-Convex, scheduled")
# #bar([1,2,3,6,7,9,10,11,14,15], [z1n, z1, z1c, z2, z2c, z3n, z3, z3c,z4,z4c] ./ 1000, color = "black", width = 0.15, label = "Obj. with balancing")
# bar([1,27], [0.0, 0.0], color = [255,255,255] ./ 255, width = 1.0)
# bar([2], [z_dc] ./ 1000, color = "coral", width = 1.0)
# bar([4,16], [z1ns, z3ns] ./ 1000, color = [127,205,187] ./ 255, width = 1.0, label = "SOC, sched.")
# bar([5,17], [z1ns, z3ns] ./ 1000, color = [29,145,192] ./ 255, width = 1.0, label = "SOC, bal.")
# bar([6,11,18,23], [z1c_d, z2c_d, z3c_d, z4c_d] ./ 1000, color = [174,1,126] ./ 255, width = 1.0, label = "Chebychev, sched.")
# bar([7,12,19,24], [z1c, z2c, z3c, z4c] ./ 1000, color = [122,1,119] ./ 255, width = 1.0, label = "Chebychev, bal.")
# bar([8,13,20,25], [z1s, z2s, z3s, z4s] ./ 1000, color = [210,220,220] ./ 255, width = 1.0, label = "Non-Convex, sched.")
# bar([9,14,21,26], [z1, z2, z3, z4] ./ 1000, color = [15,15,15] ./ 255, width = 1.0, label = "Non-Convex, bal.")
#
# legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 2)
# savefig(string("plots_final//objectives_4.8454.pdf"), format = :pdf)

#------------
#------------

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.055,0.19,0.94,0.795])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$u\$")
#ylim(bottom = 0.015, top = 20)
#ylabel("\$\\chi^{+}_{u}\$")

plot(string(u_buses), χs/γs, color = "lightgreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\beta_{u}\$")
plot(string(u_buses), σ_vec, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\sigma_{u}\$")
plot(string(u_buses), χs, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\chi_{u}\$")
plot(string(u_buses), γvec, color = "cornflowerblue", lw = 1.2, ls = "dashed", label = "\$\\chi\$")

legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9)
savefig(string("plots_final//symmetric.pdf"), format = :pdf)

#-----------------------
#-----------------------

fig = figure(figsize=(6, 2.0))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.10,0.23,0.89,0.76])
grid(linewidth = 0.1, linestyle = (0, (20, 20)), color = "lightgray", axis = "y")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4,length = 2.0)

xlabel("Balancing Framework")
#ylim(bottom = 95, top = 138)
ylim(bottom = 65, top = 100)
#ylim(bottom = 65, top = 128)
#xlim(left = 1, right = n_generators)
ylabel("\$z^{*} [1000\\\$]\$")
xticks([2.0,5.0,8.5,12.0,15.5], labels = ["DC", "S-SW","A-SW","S-N2N","A-N2N"])
#yticks([130,135,140], labels = ["130","135","140"])
#yticks([120,125,130], labels = ["120","125","130"])
#yticks([110,115,120], labels = ["110","115","120"])
#yticks([105,110,115,120,125], labels = ["105","110","115","120","125"])
yticks([75,80,85], labels = ["75","80","85"])

#bar([3,7,11,15], [z1c_d, z2c_d, z3c_d, z4c_d] ./ 1000, color = "mediumblue", width = 0.6, label = "Chebychev, scheduled")
#bar([1,9], [z1ns, z3ns] ./ 1000, color = "silver", width = 0.6, label = "SOC, scheduled")
#bar([2,6,10,14], [z1s, z2s, z3s, z4s] ./ 1000, color = "coral", width = 0.6, label = "Non-Convex, scheduled")
#bar([1,2,3,6,7,9,10,11,14,15], [z1n, z1, z1c, z2, z2c, z3n, z3, z3c,z4,z4c] ./ 1000, color = "black", width = 0.15, label = "Obj. with balancing")
bar([1,17], [0.0, 0.0], color = [255,255,255] ./ 255, width = 1.0)
bar([2], [z_dc] ./ 1000, color = [220,220,220] ./ 255, width = 1.0)
bar([4,11], [z1ns, z3ns] ./ 1000, color = [127,205,187] ./ 255, width = 1.0, label = "SOC")
#bar([5,17], [z1ns, z3ns] ./ 1000, color = [29,145,192] ./ 255, width = 1.0, label = "SOC, bal.")
#bar([6,11,18,23], [z1c_d, z2c_d, z3c_d, z4c_d] ./ 1000, color = [174,1,126] ./ 255, width = 1.0, label = "Chebychev, sched.")
bar([5,8,12,15], [z1c, z2c, z3c, z4c] ./ 1000, color = [122,1,119] ./ 255, width = 1.0, label = "Chebychev")
#bar([8,13,20,25], [z1s, z2s, z3s, z4s] ./ 1000, color = [210,220,220] ./ 255, width = 1.0, label = "Non-Convex, sched.")
bar([6,9,13,16], [z1, z2, z3, z4] ./ 1000, color = [15,15,15] ./ 255, width = 1.0, label = "Non-Convex")
annotate("infeas.", [5.0 - 0.3, 68], rotation = 90, color = [122,1,119] ./ 255)
annotate("infeas.", [8.0 - 0.3, 68], rotation = 90, color = [122,1,119] ./ 255)

legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 3)
savefig(string("plots_final//objectives_38.76338.pdf"), format = :pdf)

#------------------------------
#------------------------------

1 == 1

# fig = figure(figsize=(6, 2.6))
# rc("font", family = "serif", style = "italic", size = 14)
# rc("text", usetex = true)
# rc("lines", linewidth = 1)
#
# ax = fig.add_axes([0.10,0.18,0.89,0.8])
# grid(linewidth = 0.1, linestyle = (0, (20, 20)), color = "lightgray", axis = "y")
# ax.tick_params(direction = "in", top = true, right = true, width = 1.4,length = 2.0)
#
# xlabel("Balancing Framework")
# #ylim(bottom = 95, top = 148)
# ylim(bottom = 60, top = 110)
# #xlim(left = 1, right = n_generators)
# ylabel("\$z^{*} [1000\\\$]\$")
# xticks([2.0,6.5,12.5,18.5,24.5], labels = ["DC", "S-SW","A-SW","S-N2N","A-N2N"])
# yticks([65,70,75,80,85], labels = ["65","70","75","80","85"])
#
# #bar([3,7,11,15], [z1c_d, z2c_d, z3c_d, z4c_d] ./ 1000, color = "mediumblue", width = 0.6, label = "Chebychev, scheduled")
# #bar([1,9], [z1ns, z3ns] ./ 1000, color = "silver", width = 0.6, label = "SOC, scheduled")
# #bar([2,6,10,14], [z1s, z2s, z3s, z4s] ./ 1000, color = "coral", width = 0.6, label = "Non-Convex, scheduled")
# #bar([1,2,3,6,7,9,10,11,14,15], [z1n, z1, z1c, z2, z2c, z3n, z3, z3c,z4,z4c] ./ 1000, color = "black", width = 0.15, label = "Obj. with balancing")
# bar([1,27], [0.0, 0.0], color = [255,255,255] ./ 255, width = 1.0)
# bar([2], [z_dc] ./ 1000, color = "coral", width = 1.0)
# bar([4,16], [z1ns, z3ns] ./ 1000, color = [127,205,187] ./ 255, width = 1.0, label = "SOC, sched.")
# bar([5,17], [z1ns, z3ns] ./ 1000, color = [29,145,192] ./ 255, width = 1.0, label = "SOC, bal.")
# bar([6,11,18,23], [z1c_d, z2c_d, z3c_d, z4c_d] ./ 1000, color = [174,1,126] ./ 255, width = 1.0, label = "Chebychev, sched.")
# bar([7,12,19,24], [z1c, z2c, z3c, z4c] ./ 1000, color = [122,1,119] ./ 255, width = 1.0, label = "Chebychev, bal.")
# bar([8,13,20,25], [z1s, z2s, z3s, z4s] ./ 1000, color = [210,220,220] ./ 255, width = 1.0, label = "Non-Convex, sched.")
# bar([9,14,21,26], [z1, z2, z3, z4] ./ 1000, color = [15,15,15] ./ 255, width = 1.0, label = "Non-Convex, bal.")
# annotate("infeas.", [6.0 - 0.4, 61.5], rotation = 90, color = [174,1,126] ./ 255)
# annotate("infeas.", [7.0 - 0.4, 61.5], rotation = 90, color = [122,1,119] ./ 255)
# annotate("infeas.", [11.0 - 0.4, 61.5], rotation = 90, color = [174,1,126] ./ 255)
# annotate("infeas.", [12.0 - 0.4, 61.5], rotation = 90, color = [122,1,119] ./ 255)
#
# legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 2)
# savefig(string("plots_final//objectives_38.763.pdf"), format = :pdf)



##----------------------
##----------------------

u_buses = [i for i in 1:n_farms]

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.052,0.19,0.946,0.795])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$u\$")
ylim(bottom = -1.0, top = 7.0)
#ylabel("\$\\chi^{+}_{u}\$")

#plot(string(u_buses), χ1u, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\chi_{u}\$")
plot(string(u_buses), μ .* sqrt((1-ϵ)/ϵ), color = "teal", lw = 1.0, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\mu_u \\sqrt{\\frac{1-\\epsilon}{\\epsilon}}\$")
plot(string(u_buses), σ .* sqrt((1-ϵ)/ϵ), color = "mediumblue", lw = 1.0, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\sigma_u \\sqrt{\\frac{1-\\epsilon}{\\epsilon}}\$")
plot(string(u_buses), σ_cheb, color = "navy", lw = 1.0, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\tilde{\\sigma}_u\$")

# plot(string(u_buses), χ1up_cheb, color = "lightgreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$Chebychev, up$")
# plot(string(u_buses), χ1um_cheb, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$Chebychev, down\$")
# plot(string(u_buses), χ1up, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\chi_{u}\$")
# plot(string(u_buses), χ1um, color = "cornflowerblue", lw = 1.2, ls = "dashed", label = "\$\\chi\$")

legend(loc = "upper left", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 3)
savefig(string("plots_final//chi_mu.pdf"), format = :pdf)

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.052,0.19,0.946,0.795])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$u\$")
ylim(bottom = -600, top = 1500)
#xlim(left = 0.0, right = length(u_buses))
#ylabel("\$\\chi^{+}_{u}\$")

plot(string(u_buses), χ1up_cheb, color = "mediumblue", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi^{+}_{u}, Cheb.\$")
plot(string(u_buses), χ1um_cheb, color = "navy", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi^{-}_{u}, Cheb.\$")
plot(string(u_buses), χ1uc, color = "deeppink", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi_{u}, Cheb.\$")
plot(string(u_buses), χ1u, color = [8,8,8] ./ 255, lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi_{u}, Ncvx.\$")
plot(string(u_buses), χ1up, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi^{+}_{u}, Ncvx.\$")
plot(string(u_buses), χ1um, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi^{-}_{u}, Ncvx.\$")

# plot(string(u_buses), χ1up_cheb, color = "lightgreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$Chebychev, up$")
# plot(string(u_buses), χ1um_cheb, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$Chebychev, down\$")
# plot(string(u_buses), χ1up, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\chi_{u}\$")
# plot(string(u_buses), χ1um, color = "cornflowerblue", lw = 1.2, ls = "dashed", label = "\$\\chi\$")

legend(loc = "upper left", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 3)
savefig(string("plots_final//chi_cheb_ncvx.pdf"), format = :pdf)

# The term sum_i - 2 A_i m p_i  in the objective is negative if $ e^T m \geq 0 $. Assuming 0-mean systematically overestimates
# balancing prices when more power from uRES is available then forecasted and underestimated in the opposite case.
# This is shown in the case study and Figure (objectives) especially. Solving the full formulation slightly improves on the symmetric cases' objectives and significantly in the case of asymmetric balancing.
# Observe that assuming truncated distributions yields $\lvert\mu^{+}\rvert\geq\mu$ and  $\lvert\mu^{-}\rvert\geq\mu$.
# It can also be easily observed that the standard deviations $\sigma^{-} \leq \sigma$ and $\sigma^{+} \leq \sigma$ if the parent distribution is truncated at a point $m \in \mathbb{R}$.

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
