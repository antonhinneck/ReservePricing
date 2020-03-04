using Distributions

nrml = Distributions.Gaussian(0,1)
pdf(nrml, 0)/cdf(nrml, 0)

function millsRatio(dist::T where T <: Distributions.Distribution, x::T where T <: Real)
    return cdf(dist, x)/pdf(dist, x)
end

function splitGaussian(dist::Distributions.Gaussian, center::T where T <: Real)

    sigma = dist.σ
    stdNormal = Distributions.Gaussian(0,1)

    my_pdf = pdf(stdNormal, center)
    my_cdf = cdf(stdNormal, center)
    my_revCdf = 1 - cdf(stdNormal, center)

    ## Upper truncated gaussian
    μ_u = dist.μ + dist.σ * (my_pdf / my_revCdf)
    σ_sq_u = dist.σ^2 * (1 - center * my_pdf / my_revCdf -  (my_pdf / my_revCdf)^2)

    ## Lower truncated gaussian
    μ_l = dist.μ + dist.σ * my_pdf / my_cdf
    σ_sq_l = dist.σ^2 * (1 - center * my_pdf / my_cdf - (my_pdf / my_cdf)^2)

    return [μ_l, σ_sq_l], [μ_u, σ_sq_u]
end

function splitGaussians(means::Vector{T} where T <: Real, variances::Vector{T} where T <: Real, center::T where T <: Real)

    @assert length(means) == length(variances)
    n_dists = length(means)
    lower_μ = Vector{Float64}()
    lower_σsq = Vector{Float64}()
    upper_μ = Vector{Float64}()
    upper_σsq = Vector{Float64}()

    ζ = pi / (2 * pi - 4)
    ζ = 1 / ζ

    for i in 1:n_dists
        dist = Distributions.Gaussian(means[i], variances[i])
        lower, upper = splitGaussian(dist, center)
        push!(lower_μ, lower[1])
        push!(lower_σsq, lower[2])
        push!(upper_μ, upper[1])
        push!(upper_σsq, upper[2])
    end

    lower_μ = abs(lower_μ)
    lower_σ = abs(lower_σsq)
    upper_μ = abs(upper_μ)
    upper_σ = abs(upper_σsq)

    Σm = diagm(0 => (upper_σsq) * ζ)
    Σm_rt = sqrt(Σm)
    sm_sq = sum(Σm)
    upper = [lower_μ, Σm, Σm_rt, sm_sq]

    Σp = diagm(0 => (lower_σ))
    Σp_rt = sqrt(Σp)
    sm_sq = sum(Σp)
    lower = [lower_μ, Σp, Σp_rt, sm_sq]

    return lower, upper
end
