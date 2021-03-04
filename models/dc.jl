function build_dc(generators, buses, lines, farms; output_level = 0, pdet = nothing, adet = nothing)

        ## Model
        ##------
        m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG = output_level, MSK_IPAR_INTPNT_MAX_ITERATIONS = 1000))

        #println(string("min: ",α_min," max: ",α_max))
        ## Variables
        ##----------
        @variable(m, p[1:n_generators] >= 0)
        @variable(m, f[1:n_lines])
        @variable(m, θ[1:n_buses])
        @variable(m, α[g = 1:n_generators] >= 0)

        if pdet != nothing
                @constraint(m, eq1, p .== pdet)
        end

        if adet != nothing
                @constraint(m, eq2, α .== adet)
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

        @constraint(m, cc1[i in 1:n_generators], p[i] <= generators[i].Pgmax)
        @constraint(m, cc2[i in 1:n_generators], -p[i] <= -generators[i].Pgmin)

        ## Deterministic Costs
        ##--------------------
        @variable(m, d_con >= 0)
        @variable(m, d_lin >= 0)
        @variable(m, d_quad >= 0)

        @constraint(m, d_con == sum(generators[i].pi3 for i in 1:n_generators))
        @constraint(m, d_lin == sum(p[i] * generators[i].pi2 for i in 1:n_generators))
        @constraint(m, vec(vcat(0.5, d_quad, C_rt * p)) in RotatedSecondOrderCone())
        @expression(m, det_c, d_con + d_lin + d_quad)

        ## Objective
        ##----------
        @objective(m, Min, det_c)

        return m

end
