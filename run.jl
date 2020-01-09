using Pkg
using CSV
using LinearAlgebra, Distributions
using JuMP
using Mosek, MosekTools

cd(@__DIR__)
include("code_jl/input_dcopf.jl")
include("models/dccc.jl")

datadir = "data/ieee118"
buses, lines, generators = load_network(datadir)

x_vec = [l.x for l in lines]
X = diagm(0 => x_vec)


# Create unvertain generation (wind farms)
wp = 1.25
factor_σ =  1.25 * wp

farms = []
#                 capa    mva                   var    mva  bus
#                         Base                         Base
push!(farms, Farm(70.0  / 100 * wp, factor_σ * 7.0  / 100, 3))
push!(farms, Farm(147.0 / 100 * wp, factor_σ * 14.7 / 100, 8))
push!(farms, Farm(102.0 / 100 * wp, factor_σ * 10.2 / 100, 11))
push!(farms, Farm(105.0 / 100 * wp, factor_σ * 10.5 / 100, 20))
push!(farms, Farm(113.0 / 100 * wp, factor_σ * 11.3 / 100, 24))
push!(farms, Farm(84.0  / 100 * wp, factor_σ * 8.4  / 100, 26))
push!(farms, Farm(59.0  / 100 * wp, factor_σ * 5.9  / 100, 31))
push!(farms, Farm(250.0 / 100 * wp, factor_σ * 25.0 / 100, 38))
push!(farms, Farm(118.0 / 100 * wp, factor_σ * 11.8 / 100, 43))
push!(farms, Farm(76.0  / 100 * wp, factor_σ * 7.6  / 100, 49))
push!(farms, Farm(72.0  / 100 * wp, factor_σ * 7.2  / 100, 53))

line_limits= [ 175	175	500	175	175	175	500	500	500	175	175	175	175	175	175	175	175	175	175	175	500	175	175	175	175	175	175	175	175	175	500	500	500	175	175	500	175	500	175	175	140	175	175	175	175	175	175	175	175	500	500	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	500	175	175	500	500	500	500	500	500	500	175	175	500	175	500	175	175	500	500	175	175	175	175	175	175	175	500	175	175	175	175	175	175	500	500	175	500	500	200	200	175	175	175	500	500	175	175	500	500	500	175	500	500	175	175	175	175	175	175	175	175	175	175	200	175	175	175	175	175	175	175	175	175	500	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	500	175	175	175	500	175	175	175]

thermalLimitscale = 0.8
for i in 1:length(lines)
    lines[i].s_max = 0.99*thermalLimitscale * line_limits[i]/100
end

for (i,f) in enumerate(farms)
    push!(buses[f.bus].farm_list, i)
end

n_buses = size(buses, 1)
n_generators = size(generators, 1)
n_lines = size(lines, 1)
n_farms = size(farms, 1)

slack_bus = findall(b -> b.is_slack, buses)[1]

include("models/dccc.jl")
m_dccc = build_dccc(generators, buses, lines, farms)

optimize!(m_dccc)
objective_value(m_dccc)

alphas_sys = value.(m_dccc[:α])
dp = -dual.(m_dccc[:cc1])
dm = -dual.(m_dccc[:cc2])
λ_sys  = -dual.(m_dccc[:mc])
ηp = -dual.(m_dccc[:flowlim1])
ηm = -dual.(m_dccc[:flowlim2])

congested_lines = Vector{Int64}()

for i in 1:n_lines
    if abs(ηp[i]) >= 0.001 || abs(ηm[i]) >= 0.001
        push!(congested_lines, i)
    end
end

#include("models/dccc_fixedAlpha.jl")
#m_dccc_fixedAlpha = build_dccc_fixed(generators, buses, lines, farms, value.(m_dccc[:α]))

#optimize!(m_dccc_fixedAlpha)
#objective_value(m_dccc_fixedAlpha)

#dα = -dual.(m_dccc_fixedAlpha[:fixed])

include("models/dccc_n2n.jl")
m_dccc_n2n = build_dccc_n2n(generators, buses, lines, farms)

optimize!(m_dccc_n2n)
objective_value(m_dccc_n2n)
alphas = value.(m_dccc_n2n[:α])
alphas_i = Array{Float64, 1}(undef, n_generators)
for i in 1:n_generators
    for j in 1:n_farms
        alphas_i[i] += alphas[i, j]
    end
end
alphas_i
diff = alphas_sys .- alphas_i
diff = abs.(diff)
alphas_sys

λ_n2n  = -dual.(m_dccc_n2n[:mc])

value.(m_dccc_n2n[:α])

#using DelimitedFiles
#writedlm(string(@__DIR__,"\\alphas_sys.csv"), alphas_sys,",")
#writedlm(string(@__DIR__,"\\alphas_i.csv"), alphas_i,",")

include("plots.jl")
