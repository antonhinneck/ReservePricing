using Distributions

n = Normal(0, 0.07)

nb = 12
ns = 300000
batches = zeros(nb, ns)
for i in 1:nb
    batches[i, :] = rand(n, ns)
end

samples_pos = Array{Vector{Float64}, 1}(undef, nb)
samples_neg = Array{Vector{Float64}, 1}(undef, nb)

for i in 1:nb
    samples_pos[i] = Vector{Float64}()
    samples_neg[i] = Vector{Float64}()
    for j in 1:ns
        val = batches[i, j]
        t1 = 0
        t2 = 0
        if (val < 0.0)
            t1 = 1
        elseif (val == 0.0)
            t2 = 1
        end
        sample_pos = val * (1 - t1) * (1 - t2) + val * (t2)
        sample_neg = val * (t1) * (1 - t2) + val * (t2)
        #print("[",val,", ", sample_pos,", ", sample_neg,"]  ")
        push!(samples_pos[i], sample_pos)
        push!(samples_neg[i], sample_neg)
    end
end

μ = [mean(batches[i, :]) for i in 1:nb]
σ = [sqrt(var(batches[i, :])) for i in 1:nb]
mean(μ)
mean(σ)

μm = [mean(samples_pos[i]) for i in 1:nb]
σm = [sqrt(var(samples_pos[i])) for i in 1:nb]
mean(μm)
si = mean(σm)

si^2

μp = [mean(samples_neg[i]) for i in 1:nb]
σp = [sqrt(var(samples_neg[i])) for i in 1:nb]
mean(μp)
mean(σp)

ζ2 = (2 * pi - 4) / (2 * pi)
sp = sqrt(ζ2) * 0.07\
0.07^2
2 * sp^2

std_n = Normal(0, 1.0)
trunc_n = truncated(n, 0, Inf)

ϵ = 0.01
ρ = 0.5
z = quantile(std_n, 1 - ϵ)

## Rectified
#  (μm  + μp) + zr2 * sqrt(σm^2 + σp^2) == 0 + z * sqrt(σ)
#  => zr2 = z * σ / (sqrt(2) * mean(σm)), as σm = σp
zr2 = z * mean(σ) / (sqrt(2) * mean(σm))
mean(σ) * z
mean(μm) + mean(σm) * zr2

#------------------------------------
# # mean(μm) + mean(σm) * zr == σ * z
# # zr = (σ * z - mean(μm)) / mean(σm)
# zr = (mean(σ) * z - mean(μm)) / mean(σm)
# @assert mean(μm) + mean(σm) * zr == mean(σ) * z
