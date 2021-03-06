function build_dccc_a_apx(generators, buses, lines, farms; αm_min = nothing, αm_max = nothing, αp_min = nothing, αp_max = nothing, output_level = 0)

    ## Model
    ##------
    m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG=output_level))

    if αm_min == nothing
        αm_min = zeros(n_generators)
    end

    if αm_max == nothing
        αm_max = ones(n_generators)
    end

    if αp_min == nothing
        αp_min = zeros(n_generators)
    end

    if αp_max == nothing
        αp_max = ones(n_generators)
    end

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
    @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farmids) > 0 ? sum(farms[k].forecast for k in buses[i].farmids) : 0.0)
    @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)

    @constraint(m, B * θ .== f)
    @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].u)
    @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].u)

    @constraint(m, γp, sum(αp[i] for i in 1:n_generators) == 1)
    @constraint(m, γm, sum(αm[i] for i in 1:n_generators) == 1)

    @constraint(m, cc1[i in 1:n_generators], p[i] + μms * (αm[i] - αp[i]) + za * sm * αm[i] <= generators[i].Pgmax)
    @constraint(m, cc2[i in 1:n_generators], -p[i] + μps * (αp[i] - αm[i]) + za * sp * αp[i] <= -generators[i].Pgmin)

    @variable(m, d_con >= 0)
    @variable(m, d_lin >= 0)
    @variable(m, d_quad >= 0)

    @constraint(m, d_con == sum(generators[i].pi3 for i in 1:n_generators))
    @constraint(m, d_lin == sum(p[i] * generators[i].pi2 for i in 1:n_generators))
    @constraint(m, vec(vcat(0.5, d_quad, C_rt * p)) in RotatedSecondOrderCone())
    @expression(m, det_c, d_con + d_lin + d_quad)

    ## McCormick Envelope
    ##-------------------
    @variable(m, d_bil >= 0)
    @variable(m, ψm[1:n_generators] >= 0)
    @variable(m, ψp[1:n_generators] >= 0)

    @constraint(m, aprx11[g in 1:n_generators], ψm[g] >= αm[g] * generators[g].Pgmin + αm_min[g] * p[g] - αm_min[g] * generators[g].Pgmin)
    @constraint(m, aprx12[g in 1:n_generators], ψm[g] >= αm[g] * generators[g].Pgmax + αm_max[g] * p[g] - αm_max[g] * generators[g].Pgmax)
    @constraint(m, aprx13[g in 1:n_generators], ψm[g] <= αm[g] * generators[g].Pgmax + αm_min[g] * p[g] - αm_min[g] * generators[g].Pgmax)
    @constraint(m, aprx14[g in 1:n_generators], ψm[g] <= αm[g] * generators[g].Pgmin + αm_max[g] * p[g] - αm_max[g] * generators[g].Pgmin)

    @constraint(m, aprx21[g in 1:n_generators], ψp[g] >= αp[g] * generators[g].Pgmin + αp_min[g] * p[g] - αp_min[g] * generators[g].Pgmin)
    @constraint(m, aprx22[g in 1:n_generators], ψp[g] >= αp[g] * generators[g].Pgmax + αp_max[g] * p[g] - αp_max[g] * generators[g].Pgmax)
    @constraint(m, aprx23[g in 1:n_generators], ψp[g] <= αp[g] * generators[g].Pgmax + αp_min[g] * p[g] - αp_min[g] * generators[g].Pgmax)
    @constraint(m, aprx24[g in 1:n_generators], ψp[g] <= αp[g] * generators[g].Pgmin + αp_max[g] * p[g] - αp_max[g] * generators[g].Pgmin)

    @constraint(m, bilinear_costs,  2 * sum(generators[g].pi1 * (μms * ψm[g] - μps * ψp[g]) for g in 1:n_generators) == d_bil)

    ## Balancing Cost
    ##---------------
    @variable(m, u_quad_m >= 0)
    @variable(m, u_quad_p >= 0)
    @constraint(m, vec(vcat(0.5, u_quad_m, C_rt * αm .* sm)) in RotatedSecondOrderCone())
    @constraint(m, vec(vcat(0.5, u_quad_p, C_rt * αp .* sp)) in RotatedSecondOrderCone())
    @expression(m, unc_c, u_quad_m + u_quad_p)


    ## Objective
    ##----------
    @objective(m, Min, det_c + unc_c + d_bil)

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
