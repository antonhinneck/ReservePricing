cd(@__DIR__)
include("pkgs.jl")
include("code_jl/input.jl")
include("code_jl/gradients.jl")
#include("code_jl/TruncatedGaussian.jl")

# using Distributions
# n = Normal(0, 1)
# pdf(n, 0)
# cdf(n, 0)
#
# pdf(n, 0) / cdf(n, 0)
# 0.9/0.9999
#

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
uRESs, n_ures, σ, Σ, s_sq, Σ_rt, s, μ = create_wind_farms(fc = 2, scaling_sigma = ss)

sum([u.forecast for u in uRESs])

μs = sum(μ)
𝛭 = sum(μ[f] for f in 1:length(uRESs))

ζ1 = (sqrt(pi / 2))
ζ2 = (2 * pi - 4) / (2 * pi)

truncated_dists_p = Vector{Distribution}()
truncated_dists_m = Vector{Distribution}()
for i in 1:length(uRESs)
    push!(truncated_dists_p, TruncatedNormal(0, σ[i], 0, Inf64))
    push!(truncated_dists_m, TruncatedNormal(0, σ[i], -Inf64, 0))
end

u = 4
mean(truncated_dists_p[u])
mean(truncated_dists_m[u])
sqrt(var(truncated_dists_p[u]))
sqrt(var(truncated_dists_m[u]))

σm = [sqrt(var(d)) for d in truncated_dists_m]
σp = [sqrt(var(d)) for d in truncated_dists_p]

μm = [mean(d) for d in truncated_dists_m]
μm = μm * (-1)
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
z_cheb = sqrt((1 - ϵ)/ϵ)
# σ_cheb = sqrt((1 - ϵ)/ϵ) .* σ .- μ
# σ_cheb_m = sqrt((1 - ϵ)/ϵ) .* σm .- μm
# σ_cheb_p = sqrt((1 - ϵ)/ϵ) .* σp .- μp

# sqrt(ζ2) * s
#
# z * s
# za + sm
# sum(μ)

sum(μ) + (za + sm)
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
# counter = 1
# for f in uRESs
#     print(string(round(σ_cheb[counter], digits = 4)," & "))
#     global counter += 1
# end
#
# counter = 1
# for f in uRESs
#     print(string(round(σ_cheb_m[counter], digits = 4)," & "))
#     global counter += 1
# end
#
# counter = 1
# for f in uRESs
#     print(string(round(σ_cheb_p[counter], digits = 4)," & "))
#     global counter += 1
# end

function get_payments(m; p_det = nothing)

    generation = zeros(length(buses))
    genids = [g.busidx for g in generators]

    if p_det == nothing
        ctr = 1
        for g in genids
            generation[g] = value.(m[:p])[ctr]
            ctr += 1
        end
    else
        ctr = 1
        for g in genids
            generation[g] = p_det[ctr]
            ctr += 1
        end
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
