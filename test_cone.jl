using Mosek
using MosekTools

m = Model(with_optimizer(Mosek.Optimizer))

@variable(m, α[1:n_generators, 1:n_farms] >= 0)

@variable(m, p[1:n_generators])
@expression(m, mynorm, α * Σ_rt)
@constraint(m, cone[i in 1:n_generators], vec(vcat(p[i], mynorm[i,:])) in SecondOrderCone())

@constraint(m, adq[f in 1:n_farms], sum(α[i,f] for i in 1:n_generators) == 1)

@objective(m, Max, sum(p[i] for i in 1:n_generators))

optimize!(m)
objective_value(m)
