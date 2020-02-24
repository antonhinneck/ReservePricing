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

B = X ^ (-1) * A
B_node = A' * B

d = [b.Pd for b in buses]
sum([b.Pd for b in buses])
sum([g.Pgmax for g in generators])

include("code_jl/linApprox.jl")

include("code_jl//farms.jl")
farms, n_farms, σ_sq_vec, Σ, s_sq, Σ_rt, s = create_wind_farms()
[f.σ for f in farms]
u_buses = [f.bus for f in farms]
μ_vec = [f.μ for f in farms]
p_U = μ_vec
ν = sum(μ_vec)

include("code_jl/TruncatedGaussian.jl")
lower, upper = splitGaussians(zeros(length(μ_vec)), sqrt(σ_sq_vec), 0.0)
μm = upper[1]
Σm = upper[2]
Σm_rt = upper[3]
sm_sq = upper[4]
sm = sqrt(sm_sq)
μp = lower[1]
Σp = lower[2]
Σp_rt = lower[3]
sp_sq = lower[4]
sp = sqrt(sp_sq)
sum(Σp_rt)

for (i,f) in enumerate(farms)
    push!(buses[f.bus].farmids, i)
end

## Stochastic parameters
########################

ϵ = 0.01
z = quantile(Normal(0,1), 1-ϵ)
za = quantile(Normal(0,1), 1-ϵ/2)

## Generation Costs
###################

c_vec = [g.pi1  for g in generators]
C_mat = diagm(0 => c_vec)
C_rt = sqrt(C_mat)

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

## Models
#########

case_data, generators = updateGen(0.0, 0.48)

include("models/dccc.jl")
m_dccc = build_dccc(generators, buses, lines, farms)
optimize!(m_dccc)
termination_status(m_dccc)
z1 = objective_value(m_dccc)
a_s = value.(m_dccc[:α]) #* sum(σ_vec)
λ_s = -dual.(m_dccc[:mc])
γ = dual.(m_dccc[:γ])
p_s = value.(m_dccc[:p])

include("models/dccc_n2n.jl")
m_dccc_n2n = build_dccc_n2n(generators, buses, lines, farms)
optimize!(m_dccc_n2n)
termination_status(m_dccc_n2n)
z2 = objective_value(m_dccc_n2n)
z2_d = value.(m_dccc_n2n[:det_c])
z2_u = value.(m_dccc_n2n[:unc_c])
λ_s_n2n = -dual.(m_dccc_n2n[:mc])
χ = dual.(m_dccc_n2n[:χ])

include("models/dccc_ab.jl")
m_dccc_ab = build_dccc_ab(generators, buses, lines, farms)
optimize!(m_dccc_ab)
z3 = objective_value(m_dccc_ab)
sum(value.(m_dccc_ab[:ucp]))
sum(value.(m_dccc_ab[:cp]))
z3_up = value.(m_dccc_ab[:ucp])
z3_um = value.(m_dccc_ab[:cp])
ap = value.(m_dccc_ab[:αp]) * sum(σ_vec)
am = value.(m_dccc_ab[:αm]) * sum(σ_vec)
λ_ab  = -dual.(m_dccc_ab[:mc])
γp = dual.(m_dccc_ab[:γp])
γm = dual.(m_dccc_ab[:γm])

include("models/dccc_n2n_ab.jl")
m_dccc_n2n_ab = build_dccc_n2n_ab(generators, buses, lines, farms)
optimize!(m_dccc_n2n_ab)
termination_status(m_dccc_n2n_ab)
sum(value.(m_dccc_n2n_ab[:cp]))
sum(value.(m_dccc_n2n_ab[:ecp]))
z4 = objective_value(m_dccc_n2n_ab)
χp = dual.(m_dccc_n2n_ab[:χp])
χm = dual.(m_dccc_n2n_ab[:χm])
λ_n2n_ab = -dual.(m_dccc_n2n_ab[:mc])

## SCENARIOS σ scaling
##--------------------
#=
scenarios_chi = Vector{Vector{Float64}}()
scenarios_sigma = Vector{Vector{Float64}}()
scenarios_zu = Vector{Float64}()
scenarios_z = Vector{Float64}()
scenarios_sxs = Vector{Float64}()

scalings = [i for i in range(1, 4, step = 0.5)]

for i in scalings

    global scenario_farms, nf, σ_vec, s_sq, s_rt, s, Σ_sq = create_wind_farms(wind_buses, wind_cpcty, scaling_sigma = i, scaling_cap = 1.0)

    include("models/dccc_n2n_ab.jl")
    s_m_dccc_n2n_ab = build_dccc_n2n_ab(generators, buses, lines, scenario_farms)
    optimize!(s_m_dccc_n2n_ab)

    z = objective_value(s_m_dccc_n2n_ab)

    s_zu = value.(s_m_dccc_n2n_ab[:r_uncert])
    s_χm = dual.(s_m_dccc_n2n_ab[:χm])

    s_sxs = sum(σ_vec)

    push!(scenarios_chi, s_χm)
    push!(scenarios_sigma, σ_vec)
    push!(scenarios_zu, s_zu)
    push!(scenarios_z, z)
    push!(scenarios_sxs, s_sxs)

end

include("plots_scenarios_sigma.jl")=#

## Scenarios limits scaling
##----------------------

scenarios_chiP = Vector{Vector{Float64}}()
scenarios_chiM = Vector{Vector{Float64}}()
scenarios_dP = Vector{Vector{Float64}}()
scenarios_dM = Vector{Vector{Float64}}()
scenarios_obj = Vector{Float64}()
limits = [1.0, 0.6, 0.48]

for(i, l) in enumerate(limits)

    case_data, generators = updateGen(0.0, l)

    include("models/dccc_n2n_ab.jl")
    sp_m_dccc_n2n_ab = build_dccc_n2n_ab(generators, buses, lines, farms)
    optimize!(sp_m_dccc_n2n_ab)

    zp = objective_value(sp_m_dccc_n2n_ab)

    sp_χp = dual.(sp_m_dccc_n2n_ab[:χp])
    sp_χm = dual.(sp_m_dccc_n2n_ab[:χm])
    sp_δp = dual.(sp_m_dccc_n2n_ab[:cc1])
    sp_δm = dual.(sp_m_dccc_n2n_ab[:cc2])

    push!(scenarios_chiP, sp_χp)
    push!(scenarios_chiM, sp_χm)
    push!(scenarios_dP, sp_δp)
    push!(scenarios_dM, sp_δm)
    push!(scenarios_obj, zp)

end

include("plots_scenarios_limits.jl")

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

include("save_data.jl")
