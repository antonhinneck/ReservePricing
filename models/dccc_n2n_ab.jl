function build_dccc_n2n_ab(generators, buses, lines, farms)

    ## Model
    ##------
    output_level = 1
    m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG=output_level))

    ## Variables
    ##----------
    @variable(m, p[1:n_generators] >= 0)
    @variable(m, f[1:n_lines])
    @variable(m, θ[1:n_buses])
    @variable(m, αp[1:n_generators, 1:n_farms] >= 0)
    @variable(m, αm[1:n_generators, 1:n_farms] >= 0)

    ## General Constraints
    ##--------------------
    @constraint(m, θ[slack_bus] == 0)
    @expression(m, p_by_bus[i=1:n_buses], length(buses[i].gen_list) > 0 ? sum(p[k] for k in buses[i].gen_list) : 0.0)
    @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farm_list) > 0 ? sum(farms[k].μ for k in buses[i].farm_list) : 0.0)
    @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)

    @constraint(m, B * θ .== f)
    @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].s_max)
    @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].s_max)

    @constraint(m, χp[u in 1:n_farms], sum(αp[i, u] for i in 1:n_generators) == 1)
    @constraint(m, χm[u in 1:n_farms], sum(αm[i, u] for i in 1:n_generators) == 1)

    #@constraint(m, a[u in 1:n_farms, i in 1:n_generators], αm[i, u] == αp[i, u])

    @variable(m, pp_uncert[1:n_generators] >= 0)
    @variable(m, pm_uncert[1:n_generators] >= 0)
    @expression(m, norm_up, αp * Σ_sq)
    @expression(m, norm_dwn, αm * Σ_sq)
    @constraint(m, uncert_gen1[i in 1:n_generators], vec(vcat(pp_uncert[i], norm_up[i, :])) in SecondOrderCone())
    @constraint(m, uncert_gen2[i in 1:n_generators], vec(vcat(pm_uncert[i], norm_dwn[i, :])) in SecondOrderCone())

    @constraint(m, cc1[i in 1:n_generators], p[i] + z * pm_uncert[i] <= generators[i].g_max)
    @constraint(m, cc2[i in 1:n_generators], -p[i] + z * pp_uncert[i] <= 0)

    ## Linear Cost
    ##------------
    @variable(m, r_lin >= 0)
    @expression(m, linear_cost, sum(p[i] * generators[i].cost for i in 1:n_generators))
    @constraint(m, r_lin == linear_cost)

    ## Quadratic Cost
    ##---------------
    @variable(m, r_uncert >= 0)
    @variable(m, r_sched >= 0)
    @constraint(m, vec(vcat(0.5, r_uncert, C_rt' * pm_uncert * s)) in RotatedSecondOrderCone())
    #@constraint(m, vec(vcat(0.5, r_uncert, C_rt' * pp_uncert, C_rt' * pm_uncert)) in RotatedSecondOrderCone())
    @constraint(m, vcat(0.5, r_sched, C_rt' * p) in RotatedSecondOrderCone())
    @expression(m, quad_cost, r_sched + r_uncert)

    ## Objective
    ##----------
    @objective(m, Min, r_lin + quad_cost)

    return m

end
