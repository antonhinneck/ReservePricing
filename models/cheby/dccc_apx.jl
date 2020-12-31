function build_dccc_apx(generators, buses, lines, farms; α_min = nothing, α_max = nothing, set_sol = false, output_level = 0)

    ## Model
    ##------
    m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG = output_level, MSK_IPAR_INTPNT_MAX_ITERATIONS = 1000))

    if α_min == nothing
        α_min = zeros(n_generators)
    end

    if !set_sol && α_max == nothing
        α_max = ones(n_generators)
    end
    #println(string("min: ",α_min," max: ",α_max))
    ## Variables
    ##----------
    @variable(m, p[1:n_generators] >= 0)
    @variable(m, f[1:n_lines])
    @variable(m, θ[1:n_buses])
    @variable(m, α_min[g] <= α[g = 1:n_generators] <= α_max[g])

    if set_sol
        @assert α_min != nothing
        for g in 1:n_generators
            JuMP.fix(α[g], α_min[g]; force = true)
        end
    end


    ## General Constraints
    ##--------------------
    @constraint(m, θ[slack_bus] == 0)
    @expression(m, p_by_bus[i=1:n_buses], length(buses[i].genids) > 0 ? sum(p[k] for k in buses[i].genids) : 0.0)
    @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farmids) > 0 ? sum(farms[k].forecast for k in buses[i].farmids) : 0.0)
    @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)

    @constraint(m, B * θ .== f)
    @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].u)
    @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].u)

    @constraint(m, γ, sum(α[i] for i in 1:n_generators) == 1)

    @constraint(m, cc1[i in 1:n_generators], p[i] + α[i] * 𝛭 + z * α[i] * s <= generators[i].Pgmax)
    @constraint(m, cc2[i in 1:n_generators], -p[i] + α[i] * 𝛭 + z * α[i] * s <= -generators[i].Pgmin)

    ## Generation Cost
    ##----------------
    @variable(m, d_con >= 0)
    @variable(m, d_lin >= 0)
    @variable(m, d_quad >= 0)
    @variable(m, d_bil >= 0)
    @variable(m, ψ[1:n_generators] >= 0)

    @constraint(m, d_con == sum(generators[i].pi3 for i in 1:n_generators))
    @constraint(m, d_lin == sum(p[i] * generators[i].pi2 for i in 1:n_generators))
    @constraint(m, vec(vcat(0.5, d_quad, C_rt * p)) in RotatedSecondOrderCone())
    @expression(m, det_c, d_con + d_lin + d_quad)

    ## McCormick Envelope
    ##-------------------
    @constraint(m, bilinear_costs, sum( 2 * 𝛭 * ψ[g] * generators[g].pi1 for g in 1:n_generators) == d_bil )
    @constraint(m, aprx1[g in 1:n_generators], ψ[g] >= α[g] * generators[g].Pgmin + α_min[g] * p[g] - α_min[g] * generators[g].Pgmin)
    @constraint(m, aprx2[g in 1:n_generators], ψ[g] >= α[g] * generators[g].Pgmax + α_max[g] * p[g] - α_max[g] * generators[g].Pgmax)
    @constraint(m, aprx3[g in 1:n_generators], ψ[g] <= α[g] * generators[g].Pgmax + α_min[g] * p[g] - α_min[g] * generators[g].Pgmax)
    @constraint(m, aprx4[g in 1:n_generators], ψ[g] <= α[g] * generators[g].Pgmin + α_max[g] * p[g] - α_max[g] * generators[g].Pgmin)

    ## Balancing Cost
    ##---------------
    @variable(m, u_quad >= 0)
    @constraint(m, vec(vcat(0.5, u_quad, C_rt * α .* s)) in RotatedSecondOrderCone())
    @expression(m, unc_c, u_quad)

    ## Objective
    ##----------
    @objective(m, Min, det_c + unc_c + d_bil)

    return m

end
