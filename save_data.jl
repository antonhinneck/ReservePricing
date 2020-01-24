## EXPORT LMPs
##------------

using JLD

save("lmps//lmps.jld", "system_wide_sym", λ, "system_wide_asym", λ_ab, "n2n_sym", λ_n2n, "n2n_asym", λ_n2n_ab)

data = load("lmps//lmps.jld")
