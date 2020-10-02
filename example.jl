using Distributions

function sq(array::Array{T} where T <: Number)
    return [el^2 for el in array]
end

ϵ = 0.01
n = Normal()
z = quantile(n, 1 - ϵ)
zm = quantile(truncated(n, 0, Inf64), 1 - (ϵ / 0.5))

σ_vec = [0.07, 0.147, 0.102, 0.105, 0.113, 0.084, 0.059, 0.25, 0.118, 0.076, 0.072]

μm = Vector{Float64}()
varm = Vector{Float64}()

for σ in σ_vec
    n = Normal(0, σ)
    nt = truncated(n, 0, Inf64)
    push!(μm, mean(nt))
    push!(varm, var(nt))
end

# System-wide:

# Symmetric
# pi + ai z s ≦ Pmax
s = sqrt(sum(sq(σ_vec)))
z * s

# Asymmetric
# pi + ν ai + ai zm sm ≦ Pmax
ν = sum(μm)
sm = sqrt(sum(varm))
zm * sm

# Assuming that one generator covers all imbalance
# it hase to cover
z * s
# in the symmetric case and
z * s + ν
# in the asymmetric case.

# The same holds for the asymmetric case.
# Hence, the asymmetric cases' generation limits are more constraining,
# iff the variance of the parent distribution σ^2 > 0.
# If it is 0, the 2 models are equivalent.

# Additional test
ζ1 = (sqrt(pi / 2))
ζ2 = (2 * pi - 4) / (2 * pi)
@assert round(sm - sqrt(ζ2) * s, digits = 12) == 0
@assert all(σ_vec .* (1 / ζ1) .- μm .== 0)
