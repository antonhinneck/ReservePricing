using Distributions

nd = Normal(0, 1.0)
d = Normal(0, 0.07)
dT = TruncatedNormal(d.μ, d.σ, 0, Inf64)

ϵ = 0.01
ρ = 0.5
quantile(d, 0.999)
z = quantile(nd, 0.999)
z * 0.07

nb = 10
ns = 100000
batches = zeros(nb, ns)
for i in 1:nb
    batches[i, :] = rand(d, ns)
end

samples_pos = Array{Vector{Float64}, 1}(undef, nb)
samples_neg = Array{Vector{Float64}, 1}(undef, nb)

for i in 1:nb
    samples_pos[i] = Vector{Float64}()
    samples_neg[i] = Vector{Float64}()
    for j in 1:ns
        val = batches[i, j]
        if val >= 0
            push!(samples_pos[i], val)
        elseif val <= 0
            push!(samples_neg[i], val)
        end
    end
end

μ = [mean(batches[i, :]) for i in 1:nb]
σ = [sqrt(var(batches[i, :])) for i in 1:nb]
mean(μ)
mean(σ)

μm = [mean(samples_pos[i]) for i in 1:nb]
σm = [sqrt(var(samples_pos[i])) for i in 1:nb]
mean(μm)
mean(σm)

μp = [mean(samples_neg[i]) for i in 1:nb]
σp = [sqrt(var(samples_neg[i])) for i in 1:nb]
mean(μm)
mean(σm)
