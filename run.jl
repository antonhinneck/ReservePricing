using Pkg
using CSV
using LinearAlgebra, Distributions
using JuMP
using Mosek, MosekTools

cd(@__DIR__)
include("code_jl/input_dcopf.jl")
include("code_jl/utils.jl")
include("models/dccc.jl")

datadir = "data/ieee118"
buses, lines, generators = load_network(datadir)

x_vec = [l.x for l in lines]
X = diagm(0 => x_vec)

wind_buses = [3,8,11,20,24,26,31,38,43,49,53]
wind_cpcty = [70.0, 147.0, 102.0, 105.0, 113.0, 84.0, 59.0, 250.0, 118.0, 76.0, 72.0]

function create_wind_farms(buses::Vector{Int64}, capacity::Vector{Float64}, scaling = 1)

    @assert length(buses) == length(capacity)
    farms = Vector{Farm}()
    nf = length(buses)

    for i in 1:nf
        push!(farms,  Farm(capacity[i] / 100, scaling * capacity[i] / 10 / 100, buses[i]))
    end

    σ_vec = [i.σ for i in farms]
    s_sq = diagm(0 => (σ_vec.^2))
    s_rt = s_sq^(1/2)
    s = sum(s_rt)
    Σ = diagm(0 => (σ_vec))
    Σ_sq = sqrt(Σ)

    return farms, nf, σ_vec, s_sq, s_rt, s, Σ_sq
end

farms, nf, σ_vec, s_sq, s_rt, s, Σ_sq = create_wind_farms(wind_buses, wind_cpcty)

line_limits= [ 175	175	500	175	175	175	500	500	500	175	175	175	175	175	175	175	175	175	175	175	500	175	175	175	175	175	175	175	175	175	500	500	500	175	175	500	175	500	175	175	140	175	175	175	175	175	175	175	175	500	500	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	500	175	175	500	500	500	500	500	500	500	175	175	500	175	500	175	175	500	500	175	175	175	175	175	175	175	500	175	175	175	175	175	175	500	500	175	500	500	200	200	175	175	175	500	500	175	175	500	500	500	175	500	500	175	175	175	175	175	175	175	175	175	175	200	175	175	175	175	175	175	175	175	175	500	175	175	175	175	175	175	175	175	175	175	175	175	175	175	175	500	175	175	175	500	175	175	175]

thermalLimitscale = 0.8
for i in 1:length(lines)
    lines[i].s_max = 0.99 * thermalLimitscale * line_limits[i] / 100
end

for (i,f) in enumerate(farms)
    push!(buses[f.bus].farm_list, i)
end

n_buses = size(buses, 1)
n_generators = size(generators, 1)
n_lines = size(lines, 1)
n_farms = size(farms, 1)

slack_bus = findall(b -> b.is_slack, buses)[1]

# Stochastic parameters
#----------------------

ϵ = 0.01
z = quantile(Normal(0,1), 1-ϵ)
p_U = [f.μ for f in farms]
σ_vec = [f.σ for f in farms]

Σ = diagm(0 => (σ_vec))
Σ_sq = sqrt(Σ)

# Define B marix
# Code (based on how .index value is determined) does not work for an arbitrary pglib dataset.
# Sometimes, indexes are unique, but not in 1 ... n_lines or 1 ... n_buses.

A = zeros(n_lines, n_buses)
for i in 1:n_buses
    if size(buses[i].start_line, 1) > 0
        for l in buses[i].start_line
            A[lines[l].index, i] = 1
        end
    end
    if size(buses[i].end_line, 1) > 0
        for l in buses[i].end_line
            A[lines[l].index, i] = -1
        end
    end
end

x_vec = [l.x for l in lines]
X = diagm(0 => x_vec)

B = X ^ (-1) * A
B_node = A' * B

d = [b.d_P for b in buses]

## Generation costs
#------------------
c = [g.cost/100 for g in generators]
c_vec = [0.1 * g.cost for g in generators]
C_mat = diagm(0 => c_vec)
C_rt = sqrt(C_mat)

## MODELS, SYMMETRIC
##------------------

    ## SYMMETRIC SYSTEM-WIDE
    ##----------------------

    include("models/dccc.jl")
    m_dccc = build_dccc(generators, buses, lines, farms)

    optimize!(m_dccc)
    z1 = objective_value(m_dccc)

    z1_u = value.(m_dccc[:r_uncert])
    z1_q = value.(m_dccc[:r_sched])
    z1_l = value.(m_dccc[:r_lin])

    a_s = value.(m_dccc[:α]) * sum(σ_vec)
    λ  = -dual.(m_dccc[:mc])
    γ = -dual.(m_dccc[:γ])

    ## SYMMETRIC N2N
    ##--------------

    include("models/dccc_n2n.jl")
    m_dccc_n2n = build_dccc_n2n(generators, buses, lines, farms)

    optimize!(m_dccc_n2n)
    z2 = objective_value(m_dccc_n2n)
    termination_status(m_dccc_n2n)

    z2_u = value.(m_dccc_n2n[:r_uncert])
    z2_q = value.(m_dccc_n2n[:r_sched])
    z2_l = value.(m_dccc_n2n[:r_lin])
    sum(value.(m_dccc_n2n[:p_uncert]))
    a_n2n = value.(m_dccc_n2n[:α]) * σ_vec
    sum(value.(m_dccc_n2n[:α]))
    λ_n2n  = -dual.(m_dccc_n2n[:mc])
    χ = dual.(m_dccc_n2n[:χ])

## MODELS, ASYMMETRIC
##-------------------

    ## ASYMMETRIC SYSTEM-WIDE
    ##-----------------------

    include("models/dccc_ab.jl")
    m_dccc_ab = build_dccc_ab(generators, buses, lines, farms)
    optimize!(m_dccc_ab)
    z3 = objective_value(m_dccc_ab)

    z3_u = value.(m_dccc_ab[:r_uncert])
    z3_q = value.(m_dccc_ab[:r_sched])
    z3_l = value.(m_dccc_ab[:r_lin])
    ap = value.(m_dccc_ab[:αp]) * sum(σ_vec)
    am = value.(m_dccc_ab[:αm]) * sum(σ_vec)
    λ_ab  = -dual.(m_dccc_ab[:mc])

    γp = -dual.(m_dccc_ab[:γp])
    γm = -dual.(m_dccc_ab[:γm])

    diff = ap .- am

    ## ASYMMETRIC N2N
    ##---------------

    include("models/dccc_n2n_ab.jl")
    m_dccc_n2n_ab = build_dccc_n2n_ab(generators, buses, lines, farms)
    optimize!(m_dccc_n2n_ab)
    z4 = objective_value(m_dccc_n2n_ab)
    termination_status(m_dccc_n2n_ab)

    z4_u = value.(m_dccc_n2n_ab[:r_uncert])
    z4_q = value.(m_dccc_n2n_ab[:r_sched])
    z4_l = value.(m_dccc_n2n_ab[:r_lin])
    λ_n2n_ab  = -dual.(m_dccc_n2n_ab[:mc])
    χp = dual.(m_dccc_n2n_ab[:χp])
    χm = dual.(m_dccc_n2n_ab[:χm])
    ap_n2n = value.(m_dccc_n2n_ab[:αp]) * σ_vec
    am_n2n = value.(m_dccc_n2n_ab[:αm]) * σ_vec

    sum(value.(m_dccc_n2n_ab[:pp_uncert]))
    sum(value.(m_dccc_n2n_ab[:pm_uncert]))

## SCENARIOS σ scaling
##--------------------

scenarios_chi = Vector{Vector{Float64}}()
scenarios_sigma = Vector{Vector{Float64}}()
scenarios_zu = Vector{Float64}()
scenarios_z = Vector{Float64}()
scenarios_sxs = Vector{Float64}()

scalings = [i for i in range(1, 4, step = 0.5)]

for i in scalings

    global scenario_farms, nf, σ_vec, s_sq, s_rt, s, Σ_sq = create_wind_farms(wind_buses, wind_cpcty, i)

    include("models/dccc_n2n_ab.jl")
    s_m_dccc_n2n_ab = build_dccc_n2n_ab(generators, buses, lines, scenario_farms)
    optimize!(s_m_dccc_n2n_ab)

    z = objective_value(s_m_dccc_n2n_ab)

    s_zu = value.(s_m_dccc_n2n_ab[:r_uncert])
    s_χm = dual.(s_m_dccc_n2n_ab[:χm])

    σ_vec = [i.σ for i in scenario_farms]
    s_sxs = sum(σ_vec)

    push!(scenarios_chi, s_χm)
    push!(scenarios_sigma, σ_vec)
    push!(scenarios_zu, s_zu)
    push!(scenarios_z, z)
    push!(scenarios_sxs, s_sxs)

end

include("plots_scenarios.jl")

## EXPORT
##-------

headings1 = ["Model", "\$z^{*}\$", "\$z^{*}_{q}\$", "\$z^{*}_{u}\$", "\$z^{*}_{l}\$"]
headings2 = ["","","","", ""]
types = [Int64, Float64, Float64, Float64, Float64]
body = hcat([1, 2, 3, 4], [z1, z2, z3, z4], [z1_q, z2_q, z3_q, z4_q], [z1_u, z2_u, z3_u, z4_u], [z1_l, z2_l, z3_l, z4_l])

TexTable("texTables//objective.txt", headings1, headings2, body, types, 2)

gens = Vector{Int64}()

for i in 1:n_generators
    push!(gens, generators[i].bus_idx)
end

l1 = Vector{Float64}()
l2 = Vector{Float64}()
l3 = Vector{Float64}()
l4 = Vector{Float64}()
node_idx = Vector{Int64}()

for i in 1:n_farms
    push!(node_idx, farms[i].bus)
    push!(l1, λ[farms[i].bus] / 100)
    push!(l2, λ_n2n[farms[i].bus] / 100)
    push!(l3, λ_ab[farms[i].bus] / 100)
    push!(l4, λ_n2n_ab[farms[i].bus] / 100)
end

headings1 = ["i", "sym", "asym", "sym_n2n"]
headings2 = ["","","",""]
types = [Int64, Float64, Float64, Float64]
body = hcat(node_idx, l1, l3, l2)

TexTable("texTables//prices.txt", headings1, headings2, body, types)

headings1 = ["","Model", "sym", "asym", "asym", "sym", "asym", "asym"]
headings2 = ["\$i\$", "\$c_{i}\$", "\$\\alpha_{i}\$", "\$\\alpha^{-}_{i}\$", "\$\\alpha^{+}_{i}\$", "\$e^{T}A_{i}\$", "\$e^{T}A^{-}_{i}\$", "\$e^{T}A^{+}_{i}\$"]
types = [Int, Float64, Float64, Float64, Float64, Float64, Float64, Float64]
body = hcat(gens, c, a_s, ap, am, a_n2n, ap_n2n, am_n2n)

TexTable("texTables//alphas.txt", headings1, headings2, body, types)

s = 0
ures_cap = 0
gen_cap = 0
for i in 1:n_farms
    global ures_cap += farms[i].μ
    global s+= farms[i].σ
end
for i in 1:n_generators
    global gen_cap += generators[i].g_max
end

headings1 = ["scenario", "\$z^{*}_{\\epsilon_{g}}\$", "\$s^{2}\$", "\$\\sum_{u}\\mu_{u}\$", "\$\\sum_{i}\\bar{P}_{i}\$", "\$\\% uRES\$"]
headings2 = ["","","","", "", ""]
types = [Int64, Float64, Float64, Float64, Float64, Float64]
body = hcat([1], [z], [s^2], [ures_cap], [gen_cap], [ures_cap / (ures_cap + gen_cap)])

TexTable("texTables//system.txt", headings1, headings2, body, types, 2)


#using DelimitedFiles
#writedlm(string(@__DIR__,"\\alphas_sys.csv"), alphas_sys,",")
#writedlm(string(@__DIR__,"\\alphas_i.csv"), alphas_i,",")

include("plots.jl")
