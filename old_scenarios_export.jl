
#=
ζ2 = (2 * pi - 4) / (2 * pi)
am_n2n_ab[1,:]' * μm
ζ1 * am_n2n_ab[1,:]' * σ_vec
(pp_u[1] + pm_u[1]) / 2

pp_u[1]
pm_u[1]

norm(am_n2n_ab[1,:]' * Σ_rt) * sqrt(ζ2)
norm(am_n2n_ab[1,:]' * Σm_rt)

vecArr = Vector{Float64}()
for i in 1:length(χp)
    push!(vecArr, χp[i] + χm[i])
end
sum(vecArr) / ζ3
vcat(1, [1,2,3])
a = 0
a = 0=#

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

include("plots_scenarios_sigma.jl")

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

include("plots_scenarios_limits.jl")=#

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
