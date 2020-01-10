function build_dccc_n2n(generators, buses, lines, farms)

    # Model
    output_level = 1
    m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG=output_level))

    # Variables
    #----------
    @variable(m, p[1:n_generators] >= 0)
    @variable(m, f[1:n_lines])
    @variable(m, θ[1:n_buses])
    @variable(m, α[1:n_generators, 1:n_farms] >= 0)

    # Constraints
    #------------
    @constraint(m, θ[slack_bus] == 0)
    @expression(m, p_by_bus[i=1:n_buses], length(buses[i].gen_list) > 0 ? sum(p[k] for k in buses[i].gen_list) : 0.0)
    @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farm_list) > 0 ? sum(farms[k].μ for k in buses[i].farm_list) : 0.0)
    @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)
    #.+ pu_by_bus

    @constraint(m, B * θ .== f)
    @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].s_max)
    @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].s_max)

    @constraint(m, reserve[u in 1:n_farms], sum(α[i, u] for i in 1:n_generators) == 1)

    @variable(m, p_uncert[1:n_generators] >= 0)
    @constraint(m, conicgen[i in 1:n_generators], vec(vcat(p_uncert[i], α[i, :]' * Σ_sq)) in SecondOrderCone())

    @constraint(m, cc1[i in 1:n_generators], p[i] + z * p_uncert[i] <= generators[i].g_max)
    @constraint(m, cc2[i in 1:n_generators], -p[i] + z * p_uncert[i] <= 0)

    # Linear Cost
    @expression(m, linear_cost, sum(p[i] * generators[i].cost for i in 1:n_generators))

    # Uncertain Cost
    @variable(m, r_uncert >= 0)
    @variable(m, r_sched >= 0)
    #@constraint(m, vec(vcat(r_uncert, 0.5, C_rt * α * s * ones(n_farms))) in RotatedSecondOrderCone())
    #@constraint(m, vec(vcat(r_uncert, 0.5, c_vec' * α * Σ)) in RotatedSecondOrderCone())
    @constraint(m, cost, vec(vcat(r_uncert, 0.5, C_rt * α * Σ_sq)) in RotatedSecondOrderCone())
    @constraint(m, vcat(r_sched, 0.5, C_rt * p) in RotatedSecondOrderCone())
    @expression(m, quad_cost, r_sched + r_uncert)

    @objective(m, Min, linear_cost)

    return m

end
