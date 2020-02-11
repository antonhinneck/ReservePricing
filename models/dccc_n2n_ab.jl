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
    @expression(m, p_by_bus[i=1:n_buses], length(buses[i].genids) > 0 ? sum(p[k] for k in buses[i].genids) : 0.0)
    @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farmids) > 0 ? sum(farms[k].μ for k in buses[i].farmids) : 0.0)
    @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)

    @constraint(m, B * θ .== f)
    @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].u)
    @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].u)

    @constraint(m, χp[u in 1:n_farms], sum(αp[i, u] for i in 1:n_generators) == 1)
    @constraint(m, χm[u in 1:n_farms], sum(αm[i, u] for i in 1:n_generators) == 1)

    #@constraint(m, a[u in 1:n_farms, i in 1:n_generators], αm[i, u] == αp[i, u])

    @variable(m, pp_uncert[1:n_generators] >= 0)
    @variable(m, pm_uncert[1:n_generators] >= 0)
    @expression(m, norm_up, αp * Σ_rt)
    @expression(m, norm_dwn, αm * Σ_rt)
    @constraint(m, uncert_gen1[i in 1:n_generators], vec(vcat(pp_uncert[i], norm_up[i, :])) in SecondOrderCone())
    @constraint(m, uncert_gen2[i in 1:n_generators], vec(vcat(pm_uncert[i], norm_dwn[i, :])) in SecondOrderCone())

    @constraint(m, cc1[i in 1:n_generators], p[i] + μ_vec' * αm[i, :] + z * pm_uncert[i] <= generators[i].Pgmax)
    @constraint(m, cc2[i in 1:n_generators], -p[i] + μ_vec' * αp[i, :] + z * pp_uncert[i] <= -generators[i].Pgmin)

    @variable(m, cp[1:n_generators] >= 0)
    @variable(m, ecp[1:n_generators] >= 0)

    @constraint(m, det_approx[i in 1:n_generators, j in 1:n_coefs], cp[i] >= coefs[j][1] * p[i] + coefs[j][2])
    @constraint(m, unc_approx[i in 1:n_generators, j in 1:n_coefs], ecp[i] >= cp[i] + 0.5 * coefs[j][1] * ((αm[i, :] - αp[i, :])' * μ_vec))

    @expression(m, costs, sum(ecp[i]  for i in 1:n_generators))

    #=
    ## Linear Cost
    ##------------
    @variable(m, r1_lin >= 0)
    @expression(m, δp, αp * μ_vec)
    @expression(m, δm, αm * μ_vec)
    @expression(m, linear_cost_lin, 0.5 * sum(generators[i].pi2 * (δm[i] - δp[i]) for i in 1:n_generators))
    @constraint(m, r1_lin == linear_cost_lin)=#

    #=
    ## TODO:
    ##------
    @variable(m, r2_lin >= 0)
    @variable(m, tp0[1:n_generators] >= 0)
    @variable(m, tm0[1:n_generators] >= 0)
    @variable(m, tp1[1:n_generators] >= 0)
    @variable(m, tm1[1:n_generators] >= 0)

    @constraint(m, s1[i in 1:n_generators], 1.0 * [δp[i]; p[i]; tp0[i]]  in MOI.PowerCone(0.5))
    @constraint(m, s2[i in 1:n_generators], 1.0 * [δm[i]; p[i]; tm0[i]]  in MOI.PowerCone(0.5))
    @constraint(m, s3[i in 1:n_generators], 1.0 * [0.5, tp0[i], tp1[i]] in RotatedSecondOrderCone())
    @constraint(m, s4[i in 1:n_generators], 1.0 * [0.5, tm0[i], tm1[i]] in RotatedSecondOrderCone())
    #=
    @variable(m, r2_lin >= 0)
    @variable(m, tp[1:n_generators] >= 0)
    @variable(m, tm[1:n_generators] >= 0)
    @variable(m, tc1[1:n_generators] >= 0)
    @variable(m, tc2[1:n_generators] >= 0)
    @variable(m, tc[1:n_generators] >= 0)

    @constraint(m, sep[i in 1:n_generators], tp[i] == 0.5 * (δ[i] + p[i]))
    @constraint(m, sem[i in 1:n_generators], tm[i] == 0.5 * (δ[i] - p[i]))

    @constraint(m, sep_c1[i in 1:n_generators], vec([0.5, tc1[i], tp[i]]) in RotatedSecondOrderCone())
    @constraint(m, sep_c2[i in 1:n_generators], vec([0.5, tc2[i], tm[i]]) in RotatedSecondOrderCone())

    #@NLconstraint(m, sepm[i in 1:n_generators], tc[i] == tp[i]^2 - tm[i]^2)
    @constraint(m, sef[i in 1:n_generators], tc[i] == 0.25 * (tc1[i] - tc2[i]))=#
    #@variable(m, t[1:n_generators] >= 0)
    # δ[i] * p[i]
    #@NLexpression(m, linear_cost_lin, 0.5 * sum(generators[i].pi2 * (δ[i] * p[i]) for i in 1:n_generators))

    @constraint(m, r2_lin == 0.5 * sum(generators[i].pi2 * (tm1[i] - tp1[i]) for i in 1:n_generators))
    #@constraint(m, pAm == 0.5 * sum(generators[i].pi2 * tm1[i] for i in 1:n_generators))
    #@variable(m, r2_lin[i in 1:n_generators] >= 0)

    ## Quadratic Cost
    ##---------------
    @variable(m, r_uncert >= 0)
    @variable(m, r_sched >= 0)
    #@constraint(m, vec(vcat(0.5, r_uncert, C_rt' * pm_uncert * s)) in RotatedSecondOrderCone())
    @constraint(m, vec(vcat(0.5, r_uncert, C_rt' * pp_uncert, C_rt' * pm_uncert)) in RotatedSecondOrderCone())
    @constraint(m, vcat(0.5, r_sched, C_rt' * p) in RotatedSecondOrderCone())
    @expression(m, quad_cost, r_sched + r_uncert)

    ## Objective
    ##----------
    =#
    @objective(m, Min, costs)

    return m

end
