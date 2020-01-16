function build_dccc_ab(generators, buses, lines, farms)

    ## Model
    ##------
    output_level = 1
    m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG=output_level))


    ## Variables
    ##----------
    @variable(m, p[1:n_generators] >= 0)
    @variable(m, f[1:n_lines])
    @variable(m, θ[1:n_buses])
    @variable(m, αp[1:n_generators] >= 0)
    @variable(m, αm[1:n_generators] >= 0)

    ## General Constraints
    ##--------------------
    @constraint(m, θ[slack_bus] == 0)
    @expression(m, p_by_bus[i=1:n_buses], length(buses[i].gen_list) > 0 ? sum(p[k] for k in buses[i].gen_list) : 0.0)
    @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farm_list) > 0 ? sum(farms[k].μ for k in buses[i].farm_list) : 0.0)
    @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)

    @constraint(m, B * θ .== f)
    @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].s_max)
    @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].s_max)

    @constraint(m, γp, sum(αp[i] for i in 1:n_generators) == 1)
    @constraint(m, γm, sum(αm[i] for i in 1:n_generators) == 1)

    @constraint(m, cc1[i in 1:n_generators], p[i] + z * αm[i] * s <= generators[i].g_max)
    @constraint(m, cc2[i in 1:n_generators], -p[i] + z * αp[i] * s <= 0)

    ## Linear Cost
    ##------------
    @expression(m, linear_cost, sum(p[i] * generators[i].cost for i in 1:n_generators))

    ## Quadratic Cost
    ##---------------
    @variable(m, r_uncert >= 0)
    @variable(m, r_sched >= 0)
    @constraint(m, vec(vcat(r_uncert, 0.5, C_rt * αp * s, C_rt  * αm * s)) in RotatedSecondOrderCone())
    @constraint(m, vcat(r_sched, 0.5, C_rt * p) in RotatedSecondOrderCone())
    @expression(m, quad_cost, r_sched + r_uncert)

    ## Objective
    ##----------
    @objective(m, Min, linear_cost + quad_cost)

    return m

end
