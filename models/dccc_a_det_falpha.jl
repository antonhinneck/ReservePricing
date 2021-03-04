function build_dccc_a_det_falpha(generators, buses, lines, farms; output_level = 0, amdet = nothing, pdet = nothing, apdet = nothing, cheb = false)


        ## Model
        ##------
        m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG=output_level))

        ## Variables
        ##----------
        @variable(m, p[1:n_generators] >= 0)
        @variable(m, f[1:n_lines])
        @variable(m, θ[1:n_buses])
        # @variable(m, αp[1:n_generators] >= 0)
        # @variable(m, αm[1:n_generators] >= 0)

        # if amdet != nothing
        #         println("amdet")
        #         @constraint(m, eq21, αm .== amdet)
        # end
        #
        # if apdet != nothing
        #         println("apdet")
        #         @constraint(m, eq22, αp .== apdet)
        # end

        ## General Constraints
        ##--------------------
        @constraint(m, θ[slack_bus] == 0)
        @expression(m, p_by_bus[i=1:n_buses], length(buses[i].genids) > 0 ? sum(p[k] for k in buses[i].genids) : 0.0)
        @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farmids) > 0 ? sum(farms[k].forecast for k in buses[i].farmids) : 0.0)
        @constraint(m, mc, B_node * θ .== p_by_bus .+ pu_by_bus .- d)

        @constraint(m, B * θ .== f)
        @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].u)
        @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].u)

        # @constraint(m, γp, sum(αp[i] for i in 1:n_generators) == 1)
        # @constraint(m, γm, sum(αm[i] for i in 1:n_generators) == 1)

        # if pdet != nothing
        if !cheb
                @constraint(m, cc1[i in 1:n_generators], p[i] + μms * amdet[i] - μps * apdet[i] + za * (sm * amdet[i] + sp * apdet[i]) <= generators[i].Pgmax)
                @constraint(m, cc2[i in 1:n_generators], -p[i] + μps * apdet[i] - μms * amdet[i] + za * (sm * amdet[i] + sp * apdet[i]) <= -generators[i].Pgmin)
        else
                @constraint(m, cc1[i in 1:n_generators], p[i] + μms * amdet[i] - μps * apdet[i] + z_cheb * (sm * amdet[i] + sp * apdet[i]) <= generators[i].Pgmax)
                @constraint(m, cc2[i in 1:n_generators], -p[i] + μps * apdet[i] - μms * amdet[i] + z_cheb * (sm * amdet[i] + sp * apdet[i]) <= -generators[i].Pgmin)
        end
        # elseif amdet != nothing && apdet != nothing
        #         @constraint(m, cc1[i in 1:n_generators], p[i] + μms * (αm[i]) + za * sm * αm[i] <= generators[i].Pgmax)
        #         @constraint(m, cc2[i in 1:n_generators], -p[i] + μps * (αp[i]) + za * sp * αp[i] <= -generators[i].Pgmin)
        # else
        #         throw(Exception)
        # end

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
        # @variable(m, ψm[1:n_generators] >= 0)
        # @variable(m, ψp[1:n_generators] >= 0)
        #
        # @constraint(m, aprx11[g in 1:n_generators], ψm[g] >= αm[g] * generators[g].Pgmin + αm_min[g] * p[g] - αm_min[g] * generators[g].Pgmin)
        # @constraint(m, aprx12[g in 1:n_generators], ψm[g] >= αm[g] * generators[g].Pgmax + αm_max[g] * p[g] - αm_max[g] * generators[g].Pgmax)
        # @constraint(m, aprx13[g in 1:n_generators], ψm[g] <= αm[g] * generators[g].Pgmax + αm_min[g] * p[g] - αm_min[g] * generators[g].Pgmax)
        # @constraint(m, aprx14[g in 1:n_generators], ψm[g] <= αm[g] * generators[g].Pgmin + αm_max[g] * p[g] - αm_max[g] * generators[g].Pgmin)
        #
        # @constraint(m, aprx21[g in 1:n_generators], ψp[g] >= αp[g] * generators[g].Pgmin + αp_min[g] * p[g] - αp_min[g] * generators[g].Pgmin)
        # @constraint(m, aprx22[g in 1:n_generators], ψp[g] >= αp[g] * generators[g].Pgmax + αp_max[g] * p[g] - αp_max[g] * generators[g].Pgmax)
        # @constraint(m, aprx23[g in 1:n_generators], ψp[g] <= αp[g] * generators[g].Pgmax + αp_min[g] * p[g] - αp_min[g] * generators[g].Pgmax)
        # @constraint(m, aprx24[g in 1:n_generators], ψp[g] <= αp[g] * generators[g].Pgmin + αp_max[g] * p[g] - αp_max[g] * generators[g].Pgmin)

        @constraint(m, bilinear_costs,  2 * sum(generators[g].pi1 * p[g] * (μms * amdet[g] - μps * apdet[g]) for g in 1:n_generators) == d_bil)

        # if pdet != nothing && (amdet == nothing && apdet == nothing)
        #         println("objective --- 123")
        #         @constraint(m, bilinear_costs,  2 * sum(generators[g].pi1 * pdet[g] * (-μms * αm[g] + μps * αp[g]) for g in 1:n_generators) == d_bil)
        # elseif amdet != nothing && apdet != nothing
        #         @constraint(m, bilinear_costs,  2 * sum(generators[g].pi1 * p[g] * (-μms * amdet[g] + μps * apdet[g]) for g in 1:n_generators) == d_bil)
        # else
        #         throw(Exception)
        # end

        ## Balancing Cost
        ##---------------
        @variable(m, u_quads_m >= 0)
        @variable(m, u_quads_p >= 0)
        @constraint(m, vec(vcat(0.5, u_quads_m, C_rt * amdet .* sm)) in RotatedSecondOrderCone())
        @constraint(m, vec(vcat(0.5, u_quads_p, C_rt * apdet .* sp)) in RotatedSecondOrderCone())

        @variable(m, u_quadm_m >= 0)
        @variable(m, u_quadm_p >= 0)
        @constraint(m, vec(vcat(0.5, u_quadm_m, C_rt * amdet .* μms)) in RotatedSecondOrderCone())
        @constraint(m, vec(vcat(0.5, u_quadm_p, C_rt * apdet .* μps)) in RotatedSecondOrderCone())
        @expression(m, unc_c, u_quadm_m + u_quadm_p + u_quads_m + u_quads_p)

        ## Objective
        ##----------
        @objective(m, Min, det_c + unc_c)

        return m

end
