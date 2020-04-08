using Distributions

ζ1 = sqrt(2/pi)
ζ2 = (2 * pi - 4) / (2 * pi)

farms[1].σ

d = Normal(0, 0.07)
dT = TruncatedNormal(d.μ, d.σ, 0, Inf)

mean(dT)
std(dT)

0.07 * ζ1
0.07 * sqrt(ζ2)

[sqrt(ζ2) * f.σ for f in farms]

for i in 1:n_farms
    print(string(round([sqrt(ζ2) * f.σ for f in farms][i], digits = 3)," & "))
end

sum(μm)
sum(μm) * 0.5

tn = TruncatedNormal(0,1,0,Inf64)
za2 = quantile(tn, 1 - ϵ / 0.5)

a = farms[1].σ * z
a1 = μm[1] + sqrt(Σm[1,1]) * za
a2 = μp[1] + sqrt(Σp[1,1]) * za
pi/4
a1/a
a/a2

s * z
(sm + sum(μm)) * za / 2

d = Normal(0,1)
dt = TruncatedNormal(d.μ,d.σ,0,Inf64)

z = quantile(d, 1 - ϵ)
zt = quantile(dt, 1 - ϵ)

mean(dt) + sqrt(var(dt)) * zt
mean(d) + sqrt(var(d)) * z

ζ1 = sqrt(2/pi)
ζ2 = (2 * pi - 4) / (2 * pi)
ζ3 = (1/z) * (ζ1 * sum(σ_vec) / s + zt * sqrt(ζ2))

Σm[1,1]
sqrt(Σm[1,1])

Σ[1,1]
sqrt(Σ[1,1])

sig = sqrt(Σ[1,1])
sig / sqrt(ζ2)
μm[1]
sig / ζ1

sig * ζ3
(μm[1] + sqrt(Σm[1,1])) / ζ3

s * z
ζ3 = (1/z) * (ζ1 * sum(σ_vec) / s + zt * sqrt(ζ2))
a = 1/ζ3
b = sum(σ_vec.^2)

(sum(μm) + sm * zt) / ζ3

1/sqrt(ζ2)
1/(ζ1) * sum(σ_vec)
sum(μm)

Ai = [0.0 for i in 1:11]
Ai[1] = 0.5
Ai[2] = 0.5

ζ2 * Ai' * Σ * Ai
Ai' * Σp * Ai
sqrt(ζ2) * Ai' * Σ_rt * Ai
Ai' * Σp_rt * Ai
