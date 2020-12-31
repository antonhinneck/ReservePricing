function build_dccc_n2n(generators, buses, lines, farms; output_level = 0)

    ## Model
    ##------
    m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG=output_level))

    ## Variables
    ##----------
    @variable(m, p[1:n_generators] >= 0)
    @variable(m, f[1:n_lines])
    @variable(m, θ[1:n_buses])
    @variable(m, α[1:n_generators, 1:n_farms] >= 0)

    ## General Constraints
    ##--------------------
    @constraint(m, θ[slack_bus] == 0)
    @expression(m, p_by_bus[i=1:n_buses], length(buses[i].genids) > 0 ? sum(p[k] for k in buses[i].genids) : 0.0)
    @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farmids) > 0 ? sum(farms[k].forecast for k in buses[i].farmids) : 0.0)
    @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)

    @constraint(m, B * θ .== f)
    @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].u)
    @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].u)

    @constraint(m, χ[u in 1:n_farms], sum(α[i, u] for i in 1:n_generators) == 1)

    @variable(m, p_uncert[1:n_generators] >= 0)
    @constraint(m, cc1[i in 1:n_generators], p[i] + z * p_uncert[i] <= generators[i].Pgmax)
    @constraint(m, cc2[i in 1:n_generators], -p[i] + z * p_uncert[i] <= -generators[i].Pgmin)

    ## Deterministic Costs
    ##--------------------
    @variable(m, d_con >= 0)
    @variable(m, d_lin >= 0)
    @variable(m, d_quad >= 0)

    @constraint(m, d_con == sum(generators[i].pi3 for i in 1:n_generators))
    @constraint(m, d_lin == sum(p[i] * generators[i].pi2 for i in 1:n_generators))
    @constraint(m, vec(vcat(0.5, d_quad, C_rt * p)) in RotatedSecondOrderCone())
    @expression(m, det_c, d_con + d_lin + d_quad)

    ## Uncertain Costs
    ##----------------
    @variable(m, u_quad >= 0)
    @expression(m, norm, α * Σ_rt)
    @constraint(m, uncert_gen[i in 1:n_generators], vcat(p_uncert[i], norm[i, :]) in SecondOrderCone())

    @constraint(m, vec(vcat(0.5, u_quad, C_rt * p_uncert)) in RotatedSecondOrderCone())

    ## McCormick Envelope
    ##-------------------
    # @variable(m, u_bil >= 0)
    # #@variable(m, Ψ[1:n_generators, 1:n_farms] >= 0)
    #
    # @expression(m, μA, α_det * μ)
    # @constraint(m, bilinear_costs,  2 * sum(generators[g].pi1 * μA[g] for g in 1:n_generators) == u_bil)

    ## Quadratic
    ##----------
    @expression(m, unc_c, u_quad)

    ## Objective
    ##----------
    @objective(m, Min, unc_c + det_c)

    return m

end
