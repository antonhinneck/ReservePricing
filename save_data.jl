## EXPORT LMPs
##------------

using JLD

save("lmps//lmps.jld", "sw_sym", λ_s_c, "n2n_asym", λ_n2n_ab_c)

data = load("lmps//lmps.jld")

minimum(λ_n2n_ab_c)
