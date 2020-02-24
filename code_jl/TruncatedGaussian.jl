using Distributions

function millsRatio(dist::T where T <: Distributions.Distribution, x::T where T <: Real)
    return cdf(dist, x)/pdf(dist, x)
end

function splitGaussian(dist::Distributions.Gaussian, center::T where T <: Real)

    ## Upper truncated gaussian
    μ_u = dist.μ + dist.σ * pdf(dist, center) / (1 - cdf(dist, center))
    σ_sq_u = dist.σ^2 * (1 - center * pdf(dist, center) / (1 - cdf(dist, center)) -  (pdf(dist, center) / (1 - cdf(dist, center)))^2)

    ## Lower truncated gaussian
    μ_l = dist.μ + dist.σ * pdf(dist, center) / cdf(dist, center)
    σ_sq_l = dist.σ^2 * (1 - center * pdf(dist, center) / cdf(dist, center) -  (pdf(dist, center) / cdf(dist, center))^2)

    return [μ_l, σ_sq_l], [μ_u, σ_sq_u]
end

function splitGaussians(means::Vector{T} where T <: Real, variances::Vector{T} where T <: Real, center::T where T <: Real)

    @assert length(means) == length(variances)
    n_dists = length(means)
    lower_μ = Vector{Float64}()
    lower_σ = Vector{Float64}()
    upper_μ = Vector{Float64}()
    upper_σ = Vector{Float64}()

    for i in 1:n_dists
        dist = Distributions.Gaussian(means[i], variances[i])
        lower, upper = splitGaussian(dist, center)
        push!(lower_μ, lower[1])
        push!(lower_σ, lower[2])
        push!(upper_μ, upper[1])
        push!(upper_σ, upper[2])
    end

    lower_μ = abs(lower_μ)
    lower_σ = abs(lower_σ)
    upper_μ = abs(upper_μ)
    upper_σ = abs(upper_σ)

    Σm = diagm(0 => (upper_σ.^2))
    Σm_rt = sqrt(Σm)
    sm_sq = sum(Σm)
    upper = [lower_μ, Σm, Σm_rt, sm_sq]

    Σp = diagm(0 => (lower_σ.^2))
    Σp_rt = sqrt(Σp)
    sm_sq = sum(Σp)
    lower = [lower_μ, Σp, Σp_rt, sm_sq]

    return lower, upper
end
