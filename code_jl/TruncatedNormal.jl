dist = Distributions.Gaussian(0, 12)
quantiles = [0.25, 0.75]
truncated_means = Array{Float64, 1}(undef, length(quantiles))
cdf(normal, -1)
for q in quantiles
    μ_t =
    iters = 0
    error_limit = 10e-6
    error = 1
    y = normal.σ
    while error_limit < error && iters < 100

        μ_t = abs(cdf(normal, μ_t) - q) -
        abs(cdf(normal, μ_t) - q)

function argCDF(dist::Distributions.Gaussian, fx::Float64)



end

quantile(dist, 0.25)
cdf(dist, -8.093)
