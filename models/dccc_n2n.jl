function build_dccc_n2n(generators, buses, lines, farms)

    ## Model
    ##------
    output_level = 1
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
    @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farmids) > 0 ? sum(farms[k].μ for k in buses[i].farmids) : 0.0)
    @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)

    @constraint(m, B * θ .== f)
    @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].u)
    @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].u)

    @constraint(m, χ[u in 1:n_farms], sum(α[i, u] for i in 1:n_generators) == 1)

    @variable(m, p_uncert[1:n_generators] >= 0)
    #@variable(m, alpha_rt[1:n_generators, 1:n_farms] >= 0)

    #@constraint(m, root[i in 1:n_generators, f in 1:n_farms], 1.0 * [0.5, alpha_rt[i, f], α[i, f]]  in RotatedSecondOrderCone())

    @expression(m, norm, α * Σ) #α

    #@constraint(m, uncert_gen[i in 1:n_generators], vec(vcat(p_uncert[i], norm[i, :])) in SecondOrderCone())
    @constraint(m, uncert_gen[i in 1:n_generators], p_uncert[i] == sum(norm[i, :]))

    @constraint(m, cc1[i in 1:n_generators], p[i] + z * p_uncert[i] <= generators[i].Pgmax)
    @constraint(m, cc2[i in 1:n_generators], -p[i] + z * p_uncert[i] <= -generators[i].Pgmin)

    @variable(m, cp[1:n_generators] >= 0)

    @constraint(m, det_approx[i in 1:n_generators, j in 1:length(my_aprxs[i].coefs)], cp[i] >= my_aprxs[i].coefs[j][1] * p[i] + my_aprxs[i].coefs[j][2])

    @expression(m, costs, sum(cp[i]  for i in 1:n_generators))

    ## Objective
    ##----------

    @objective(m, Min, costs)

    ## Generation Cost
    ##----------------
    #=
    @variable(m, d_con >= 0)
    @variable(m, d_lin >= 0)
    @variable(m, d_quad >= 0)
    @constraint(m, d_con == sum(generators[i].pi3 for i in 1:n_generators))
    @constraint(m, d_lin == sum(p[i] * generators[i].pi2 for i in 1:n_generators))
    @constraint(m, vec(vcat(0.5, d_quad, C_rt * p)) in RotatedSecondOrderCone())
    @expression(m, det_c, d_con + d_lin + d_quad)

    ## Balancing Cost
    ##---------------

    @variable(m, u_quad >= 0)
    @constraint(m, vec(vcat(0.5, u_quad, C_rt * p_uncert)) in RotatedSecondOrderCone())
    @expression(m, unc_c, u_quad)

    ## Objective
    ##----------
    @objective(m, Min, unc_c + det_c)=#

    return m

end
