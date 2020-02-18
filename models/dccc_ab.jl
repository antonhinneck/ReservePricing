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
    @expression(m, p_by_bus[i=1:n_buses], length(buses[i].genids) > 0 ? sum(p[k] for k in buses[i].genids) : 0.0)
    @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farmids) > 0 ? sum(farms[k].μ for k in buses[i].farmids) : 0.0)
    @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)

    @constraint(m, B * θ .== f)
    @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].u)
    @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].u)

    @constraint(m, γp, sum(αp[i] for i in 1:n_generators) == 1)
    @constraint(m, γm, sum(αm[i] for i in 1:n_generators) == 1)

    @constraint(m, cc1[i in 1:n_generators], p[i] + z * s * αm[i] <= generators[i].Pgmax)
    @constraint(m, cc2[i in 1:n_generators], -p[i] + z * s * αp[i] <= generators[i].Pgmin)

    @variable(m, cp[1:n_generators] >= 0)
    @constraint(m, approx1[i in 1:n_generators, j in 1:n_coefs], cp[i] >= my_aprxs[i].coefs[j][1] * p[i] + my_aprxs[i].coefs[j][2])
    @expression(m, det_c, sum(cp[i] for i in 1:n_generators))

    @variable(m, ucp[1:n_generators] >= 0)
    @constraint(m, approx2[i in 1:n_generators, j in 1:n_coefs], ucp[i] >= 0.5 * ν * my_aprxs[i].coefs[j][1] * (αm[i] - αp[i]))
    @expression(m, unc_c, sum(ucp[i] for i in 1:n_generators))
    #=@variable(m, ucp[1:n_generators] >= 0)
    @constraint(m, approx2[i in 1:n_generators, j in 1:n_coefs], ucp[i] >= 0.5 * ν * coefs[j][1] * (αm[i] - αp[i]))
    @expression(m, det_c, sum(cp[i] for i in 1:n_generators))
    @expression(m, unc_c, sum(ucp[i] for i in 1:n_generators))=#

    ## Objective
    ##----------
    @objective(m, Min, det_c + unc_c)


    return m

end

#ν * αp[i] +
#=
## Generation Cost
##----------------
@variable(m, d_con >= 0)
@variable(m, d_lin >= 0)
@variable(m, d_quad >= 0)
@constraint(m, d_con == sum(generators[i].pi3 for i in 1:n_generators))
@constraint(m, d_lin == sum(p[i] * generators[i].pi2 for i in 1:n_generators))
@constraint(m, vec(vcat(0.5, d_quad, C_rt * p)) in RotatedSecondOrderCone())
@expression(m, det_c, d_con + d_lin + d_quad)

## Linear Cost
##------------
@variable(m, u_lin)
@constraint(m, u_lin == ν * sum(generators[i].pi2 * (αm[i] - αp[i]) for i in 1:n_generators))

## Bilinear Cost
##--------------

@variable(m, u_bilin)
@variable(m, tp0[1:n_generators] >= 0)
@variable(m, tm0[1:n_generators] >= 0)
@variable(m, tp1[1:n_generators] >= 0)
@variable(m, tm1[1:n_generators] >= 0)

@expression(m, δp, αp .* ν .* [generators[i].pi2 for i in 1:n_generators])
@expression(m, δm, αm .* ν .* [generators[i].pi2 for i in 1:n_generators])

@constraint(m, s1[i in 1:n_generators], 1.0 * [δp[i]; p[i]; tp0[i]]  in MOI.PowerCone(0.5))
@constraint(m, s2[i in 1:n_generators], 1.0 * [δm[i]; p[i]; tm0[i]]  in MOI.PowerCone(0.5))
@constraint(m, s3[i in 1:n_generators], 1.0 * [0.5, tp0[i], tp1[i]] in RotatedSecondOrderCone())
@constraint(m, s4[i in 1:n_generators], 1.0 * [0.5, tm0[i], tm1[i]] in RotatedSecondOrderCone())

@constraint(m, u_bilin == sum(tm1[i] - tp1[i] for i in 1:n_generators))

## Quadratic Cost
##---------------
@variable(m, u_quad_p >= 0)
@variable(m, u_quad_m >= 0)
@constraint(m, vec(vcat(0.25, u_quad_p, C_rt * αp .* s)) in RotatedSecondOrderCone())
@constraint(m, vec(vcat(0.25, u_quad_m, C_rt * αm .* s)) in RotatedSecondOrderCone())

#@constraint(m, vcat(r_sched, 0.5, C_rt * p) in RotatedSecondOrderCone())
@expression(m, unc_c, u_lin + u_quad_p + u_quad_m)=#
