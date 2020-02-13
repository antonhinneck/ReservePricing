function build_dccc_det(generators, buses, lines, farms)

    ## Model
    ##------
    output_level = 1
    m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG=output_level))

    ## Variables
    ##----------
    @variable(m, p[1:n_generators] >= 0)
    @variable(m, f[1:n_lines])
    @variable(m, θ[1:n_buses])

    ## General Constraints
    ##--------------------
    @constraint(m, θ[slack_bus] == 0)
    @expression(m, p_by_bus[i=1:n_buses], length(buses[i].genids) > 0 ? sum(p[k] for k in buses[i].genids) : 0.0)
    @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farmids) > 0 ? sum(farms[k].μ for k in buses[i].farmids) : 0.0)
    @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)

    @constraint(m, B * θ .== f)
    @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].u)
    @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].u)

    @constraint(m, cc1[i in 1:n_generators], p[i] <= generators[i].Pgmax)
    @constraint(m, cc2[i in 1:n_generators], -p[i] <= -generators[i].Pgmin)

    ## Objective
    ##----------

    @variable(m, cp[1:n_generators] >= 0)

    @constraint(m, det_approx[i in 1:n_generators, j in 1:n_coefs], cp[i] >= coefs[j][1] * p[i] + coefs[j][2])

    @expression(m, costs, sum(cp[i]  for i in 1:n_generators))

    @objective(m, Min, costs)

    return m

end
