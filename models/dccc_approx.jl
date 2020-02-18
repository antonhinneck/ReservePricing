function build_dccc_prx(generators, buses, lines, farms)

    ## Model
    ##------
    output_level = 1
    m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG=output_level))

    ## Variables
    ##----------
    @variable(m, p[1:n_generators] >= 0)
    @variable(m, f[1:n_lines])
    @variable(m, θ[1:n_buses])
    @variable(m, α[1:n_generators] >= 0)

    ## General Constraints
    ##--------------------
    @constraint(m, θ[slack_bus] == 0)
    @expression(m, p_by_bus[i=1:n_buses], length(buses[i].genids) > 0 ? sum(p[k] for k in buses[i].genids) : 0.0)
    @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farmids) > 0 ? sum(farms[k].μ for k in buses[i].farmids) : 0.0)
    @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)

    @constraint(m, B * θ .== f)
    @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].u)
    @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].u)

    @constraint(m, γ, sum(α[i] for i in 1:n_generators) == 1)

    @constraint(m, cc1[i in 1:n_generators], p[i] + z * α[i] * s <= generators[i].Pgmax)
    @constraint(m, cc2[i in 1:n_generators], -p[i] + z * α[i] * s <= -generators[i].Pgmin)

    #@variable(m, cp[1:n_generators] >= 0)
    #@constraint(m, approx[i in 1:n_generators, j in 1:n_coefs], cp[i] >= coefs[j][1] * p[i] + coefs[j][2])

    #@expression(m, det_c, sum(cp[i] for i in 1:n_generators))


    ## Generation Cost
    ##----------------

    @variable(m, cp[1:n_generators] >= 0)
    @constraint(m, approx[i in 1:n_generators, j in 1:n_coefs], cp[i] >= my_aprxs[i].coefs[j][1] * p[i] + my_aprxs[i].coefs[j][2])
    @expression(m, det_c, sum(cp[i] for i in 1:n_generators))

    ## Objective
    ##----------
    @objective(m, Min, det_c)

    return m

end
