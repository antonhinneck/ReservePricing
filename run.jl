cd(@__DIR__)
include("pkgs.jl")
include("code_jl/input.jl")

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
farms, n_farms, σ_vec, Σ, s_sq, Σ_rt, s = create_wind_farms()

u_buses = [f.bus for f in farms]
μ_vec = [f.μ for f in farms]
p_U = μ_vec
ν = sum(μ_vec)

for (i,f) in enumerate(farms)
    push!(buses[f.bus].farmids, i)
end

## Stochastic parameters
##----------------------

ϵ = 0.01
z = quantile(Normal(0,1), 1-ϵ)

## Generation costs
##-----------------

c_vec = [g.pi1  for g in generators]
C_mat = diagm(0 => c_vec)
C_rt = sqrt(C_mat)

##-----------------
## Models
##-----------------
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

case_data, generators = updateGen(0.2, 0.9)

# z = 0.01, dccc +- 0.1591
# z = 0.01, dccc_n2n +- 0.05

include("models/dccc.jl")
m_dccc = build_dccc(generators, buses, lines, farms)
optimize!(m_dccc)
z1 = objective_value(m_dccc)
termination_status(m_dccc)
z1_u = value.(m_dccc[:det_c])
z1_u = value.(m_dccc[:unc_c])
value.(m_dccc[:p])
a_s = value.(m_dccc[:α]) #* sum(σ_vec)
λ = -dual.(m_dccc[:mc])
γ = dual.(m_dccc[:γ])

include("models/dccc_n2n.jl")
m_dccc_n2n = build_dccc_n2n(generators, buses, lines, farms)
optimize!(m_dccc_n2n)
z2 = objective_value(m_dccc_n2n)
termination_status(m_dccc_n2n)
sum(dual.(m_dccc_n2n[:χ]))
dual.(m_dccc_n2n[:χ])
value.(m_dccc_n2n[:unc_c])


include("models/dccc_ab.jl")
m_dccc_ab = build_dccc_ab(generators, buses, lines, farms)
optimize!(m_dccc_ab)
dual.(m_dccc_ab[:γp])
z3 = objective_value(m_dccc_ab)
z3_up = value.(m_dccc_ab[:det_c])
z3_um = value.(m_dccc_ab[:unc_c])
ap = value.(m_dccc_ab[:αp]) * sum(σ_vec)
am = value.(m_dccc_ab[:αm]) * sum(σ_vec)
λ_ab  = -dual.(m_dccc_ab[:mc])

include("models/dccc_n2n_ab.jl")
m_dccc_n2n_ab = build_dccc_n2n_ab(generators, buses, lines, farms)
optimize!(m_dccc_n2n_ab)
z4 = objective_value(m_dccc_n2n_ab)
termination_status(m_dccc_n2n_ab)
χp = dual.(m_dccc_n2n_ab[:χp])
χm = dual.(m_dccc_n2n_ab[:χm])
#dual.(m_dccc_n2n_ab[:χm])

value.(m_dccc_n2n_ab[:cp])
value.(m_dccc_n2n_ab[:ecp])

print(unconstrained_generation)

#[generators[i].Pgmax for i in 1:n_generators]
## MODELS
##------------------

## SYMMETRIC SYSTEM-WIDE
##----------------------
#=
include("models/dccc.jl")
m_dccc = build_dccc(generators, buses, lines, farms)

optimize!(m_dccc)
z1 = objective_value(m_dccc)

z1_u = value.(m_dccc[:det_c])
#z1_q = value.(m_dccc[:unc_c])
value.(m_dccc[:p])

a_s = value.(m_dccc[:α]) #* sum(σ_vec)
λ = -dual.(m_dccc[:mc])
γ = dual.(m_dccc[:γ])

## ASYMMETRIC SYSTEM-WIDE
##-----------------------
[generators[i].pi2 for i in 1:n_generators]
[generators[i].pi1 for i in 1:n_generators]
include("models/dccc_ab.jl")
m_dccc_ab = build_dccc_ab(generators, buses, lines, farms)
optimize!(m_dccc_ab)
z3 = objective_value(m_dccc_ab)

z3_up = value.(m_dccc_ab[:det_c])
z3_um = value.(m_dccc_ab[:unc_c])
#z3_q = value.(m_dccc_ab[:r_sched])
# z3_l = value.(m_dccc_ab[:u_lin])
# z3_bl = value.(m_dccc_ab[:u_bilin])
# z3_q = value.(m_dccc_ab[:u_quad_p])
# z3_q = value.(m_dccc_ab[:u_quad_m])
ap = value.(m_dccc_ab[:αp]) * sum(σ_vec)
am = value.(m_dccc_ab[:αm]) * sum(σ_vec)
λ_ab  = -dual.(m_dccc_ab[:mc])

γp = dual.(m_dccc_ab[:γp])
γm = dual.(m_dccc_ab[:γm])
#m_dccc_ab[:δp]
diff = ap .- am
#=
@inline function evalGap(m::JuMP.Model)

    act_p = value.(m[:tp0])
    act_m = value.(m[:tm0])

    ub_p = value.(m[:δp]) .* value.(m[:p])
    ub_m = value.(m[:δm]) .* value.(m[:p])

    return act_p, ub_p, act_m, ub_m

end

evalGap(m_dccc_ab)=#

## NODE-TO-NODE
##-------------

## SYM
##----

include("models/dccc_n2n.jl")
m_dccc_n2n = build_dccc_n2n(generators, buses, lines, farms)

optimize!(m_dccc_n2n)
z2 = objective_value(m_dccc_n2n)
termination_status(m_dccc_n2n)

dual.(m_dccc_n2n[:χ])
value.(m_dccc_n2n[:p_uncert])
sum(value.(m_dccc_n2n[:α])[:,1])

#z2_u = value.(m_dccc_n2n[:r_uncert])
#z2_q = value.(m_dccc_n2n[:r_sched])
#z2_l = value.(m_dccc_n2n[:r_lin])
sum(value.(m_dccc_n2n[:p_uncert]))
a_n2n = value.(m_dccc_n2n[:α]) * σ_vec
sum(value.(m_dccc_n2n[:α]))
λ_n2n  = -dual.(m_dccc_n2n[:mc])
χ = dual.(m_dccc_n2n[:χ])
sum(χ)

## ASYM
##-----



for i in 1:n_generators
    generators[i].Pgmin = deterministic_generation[i] * 0.5
    #generators[i].Pgmax = deterministic_generation[i] * 1.4
end

include("models/dccc_n2n_ab.jl")
m_dccc_n2n_ab = build_dccc_n2n_ab(generators, buses, lines, farms)
optimize!(m_dccc_n2n_ab)
z4 = objective_value(m_dccc_n2n_ab)
termination_status(m_dccc_n2n_ab)

d1 = dual.(m_dccc_n2n_ab[:χp])
d2 = dual.(m_dccc_n2n_ab[:χm])
cp = value.(m_dccc_n2n_ab[:cp])
ecp = value.(m_dccc_n2n_ab[:ecp])

d = abs.(d1) .- abs.(d2)
c = ecp .- cp
maximum(c)

## SCENARIOS σ scaling
##--------------------

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

    #σ_vec = [i.σ for i in scenario_farms]
    s_sxs = sum(σ_vec)

    push!(scenarios_chi, s_χm)
    push!(scenarios_sigma, σ_vec)
    push!(scenarios_zu, s_zu)
    push!(scenarios_z, z)
    push!(scenarios_sxs, s_sxs)

end

include("plots_scenarios_sigma.jl")

## Scenarios pen scaling
##----------------------

scenarios_chi = Vector{Vector{Vector{Float64}}}()
scenarios_sigma = Vector{Vector{Vector{Float64}}}()
scenarios_zu =  Vector{Vector{Float64}}()
scenarios_z =  Vector{Vector{Float64}}()
scenarios_sxs = Vector{Vector{Float64}}()

sigmas = [1.0, 2.0, 3.0, 4.0]
scalings = [i for i in range(1, 2.0, step = 0.2)]

scenarios = [string(scalings[i]) for i in 1:length(scalings)]

for j in 1:length(sigmas)

    push!(scenarios_chi, Vector{Vector{Float64}}())
    push!(scenarios_sigma, Vector{Vector{Float64}}())
    push!(scenarios_zu, Vector{Float64}())
    push!(scenarios_z, Vector{Float64}())
    push!(scenarios_sxs, Vector{Float64}())

    for i in scalings

        global scenario_farms, nf, σ_vec, s_sq, s_rt, s, Σ_sq = create_wind_farms(wind_buses, wind_cpcty, scaling_sigma = sigmas[j], scaling_cap = i)

        include("models/dccc_n2n_ab.jl")
        sp_m_dccc_n2n_ab = build_dccc_n2n_ab(generators, buses, lines, scenario_farms)
        optimize!(sp_m_dccc_n2n_ab)

        zp = objective_value(sp_m_dccc_n2n_ab)

        sp_zu = value.(sp_m_dccc_n2n_ab[:r_uncert])
        sp_χm = dual.(sp_m_dccc_n2n_ab[:χm])

        σ_vec = [i.σ for i in scenario_farms]
        sp_sxs = sum(σ_vec)

        push!(scenarios_chi[j], sp_χm)
        push!(scenarios_sigma[j], σ_vec)
        push!(scenarios_zu[j], sp_zu)
        push!(scenarios_z[j], zp)
        push!(scenarios_sxs[j], sp_sxs)

    end
end

include("plots_scenarios_pen.jl")

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

include("save_data.jl")=#
