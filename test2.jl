ζ1 = (sqrt(2 / pi))
ζ2 = (2 * pi) / (2 * pi - 4)
@assert round(sm - sqrt(ζ2) * s, digits = 12) == 0
@assert all(σ_vec .* (1 / ζ1) .- μm .== 0)

σ / sqrt(ζ2)

using Gurobi
using JuMP

σ = 0.07

ζ1 = (sqrt(pi / 2))
ζ2 = (2 * pi) / (2 * pi - 4)

m = Model(Gurobi.Optimizer)

@variable(m, μm)
@variable(m, σm)

@constraint(m, σ == μm * ζ1)
@constraint(m, σ == σm * sqrt(ζ2))



optimize!(m)

value(m[:μm])
value(m[:σm])
