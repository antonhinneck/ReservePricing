function build_dccc_det(generators, buses, lines, farms, α_det::Array{T, 1} where T <: Real, p_det::Array{T, 1} where T <: Real; output_level = 0)

    ## Model
    ##------
    m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG=output_level, MSK_IPAR_INTPNT_MAX_ITERATIONS = 500))

    ## Variables
    ##----------
    @variable(m, p[1:n_generators] >= 0)
    @variable(m, f[1:n_lines])
    @variable(m, θ[1:n_buses])
    @variable(m, α[1:n_generators] >= 0)

    for g in 1:n_generators
        JuMP.fix(α[g], α_det[g]; force = true)
        #JuMP.fix(p[g], p_det[g]; force = true)
    end

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

    @constraint(m, cc1[i in 1:n_generators], p[i] + α[i] * μs + z * α[i] * s <= generators[i].Pgmax)
    @constraint(m, cc2[i in 1:n_generators], -p[i] + α[i] * μs + z * α[i] * s <= -generators[i].Pgmin)

    ## Generation Cost
    ##----------------
    @variable(m, d_con >= 0)
    @variable(m, d_lin >= 0)
    @variable(m, d_quad >= 0)
    @variable(m, d_bil >= 0)

    @constraint(m, d_con == sum(generators[i].pi3 for i in 1:n_generators))
    @constraint(m, d_lin == sum(p[i] * generators[i].pi2 for i in 1:n_generators))
    @constraint(m, vec(vcat(0.5, d_quad, C_rt * p)) in RotatedSecondOrderCone())
    @constraint(m, d_bil_c, sum(2 * generators[g].pi1 * α_det[g] * p[g] for g in 1:n_generators) == d_bil)
    @expression(m, det_c, d_con + d_lin + d_quad + d_bil)

    ## Balancing Cost
    ##---------------
    @variable(m, u_quad >= 0)
    @constraint(m, vec(vcat(0.5, u_quad, C_rt * α .* s)) in RotatedSecondOrderCone())
    @expression(m, unc_c, u_quad)

    ## Objective
    ##----------
    @objective(m, Min, det_c + unc_c)

    return m

end
