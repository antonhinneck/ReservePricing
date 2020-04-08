ζ1 = sqrt(2/pi)
ζ2 = (2*pi - 4)/(2*pi)

function compute_ζ3n2n(α)
    ζ3i = [(ζ1 * α[i,:]' * σ_vec) / (z * norm(α[i,:]' * Σ_rt)) + sqrt(ζ2) * (za/z) for i in 1:n_generators]
    return ζ3i
end

ζ3sw = (1/z) * (ζ1 * sum(σ_vec) / s + sqrt(ζ2) * za)
print(ζ3sw)
ζ3n2n = mean(compute_ζ3n2n(α))
print(minimum(compute_ζ3n2n(α)))
print(maximum(compute_ζ3n2n(α)))
print(mean(compute_ζ3n2n(α)))
