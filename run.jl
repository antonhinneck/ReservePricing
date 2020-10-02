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

include("code_jl/linApprox.jl")

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

include("plotDists.jl")

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

## Experiments
##############

## Linear Costs
###############

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
include("code_jl/linApprox.jl")

include("models/dccc_sym.jl")
m_dccc = build_dccc_sym(generators, buses, lines, farms)
optimize!(m_dccc)
termination_status(m_dccc)
z1 = objective_value(m_dccc)

include("models/dccc_sym_lin.jl")
m_dccc = build_dccc_sym_lin(generators, buses, lines, farms)
optimize!(m_dccc)
termination_status(m_dccc)
z2 = objective_value(m_dccc)

results_approx = Vector{Float64}()
my_aprxs = Vector{aprx}()

for i in 1:16

    global my_aprxs = Vector{aprx}()
    for j in 1:n_generators
        push!(my_aprxs, aprx(approx(generators, j, segments_pu = i)...))
    end

    # To approximate the quadratic cost function, linear segments are used.
    #     These are added as constraints. As these are greater than the function

    include("models/dccc_sym_lin.jl")
    m_dccc_lin = build_dccc_sym_lin(generators, buses, lines, farms)
    optimize!(m_dccc_lin)
    termination_status(m_dccc_lin)
    push!(results_approx, objective_value(m_dccc_lin))

end

include("plot_approx.jl")

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
