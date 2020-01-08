function build_dccc_fixed(generators, buses, lines, farms, alpha)

    # DC OPF with explicit theta
    # Deterministic model with wind generation

    output_level = 0
    m_det = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG=output_level))

    # Stochastic parameters
    #----------------------

    ϵ = 0.01
    z = quantile(Normal(0,1), 1-ϵ)
    p_U = [f.μ for f in farms]
    σ_vec = [f.σ for f in farms]
    s_sq = diagm(0 => (σ_vec.^2))
    s_rt = s_sq^(1/2)
    s = sqrt(sum(s_rt))

    # Define B marix
    # Code (based on how .index value is determined) does not work for an arbitrary pglib dataset.
    # Sometimes, indexes are unique, but not in 1 ... n_lines or 1 ... n_buses.

    A = zeros(n_lines, n_buses)
    for i in 1:n_buses
        if size(buses[i].start_line, 1) > 0
            for l in buses[i].start_line
                A[lines[l].index, i] = 1
            end
        end
        if size(buses[i].end_line, 1) > 0
            for l in buses[i].end_line
                A[lines[l].index, i] = -1
            end
        end
    end

    x_vec = [l.x for l in lines]
    X = diagm(0 => x_vec)

    B = X ^ (-1) * A
    B_node = A' * B

    d = [b.d_P for b in buses]

    c_vec = [g.cost * 0.1 for g in generators]
    C_mat = diagm(0 => c_vec)
    C_rt = C_mat ^ (-1/2)

    # Model
    output_level = 0
    m = Model(with_optimizer(Mosek.Optimizer,  MSK_IPAR_LOG=output_level))

    # Variables
    #----------
    @variable(m, p[1:n_generators])
    @variable(m, f[1:n_lines])
    @variable(m, θ[1:n_buses])
    @variable(m, α[1:n_generators])

    # Constraints
    #------------
    @constraint(m, θ[slack_bus] == 0)
    @expression(m, p_by_bus[i=1:n_buses], length(buses[i].gen_list) > 0 ? sum(p[k] for k in buses[i].gen_list) : 0.0)
    @expression(m, pu_by_bus[i=1:n_buses], length(buses[i].farm_list) > 0 ? sum(farms[k].μ for k in buses[i].farm_list) : 0.0)
    @constraint(m, B_node * θ .== p_by_bus .+ pu_by_bus .- d)
    @constraint(m, B * θ .== f)
    @constraint(m, flowlim1[i in 1:n_lines], f[i] <= lines[i].s_max)
    @constraint(m, flowlim2[i in 1:n_lines], -f[i] >= -lines[i].s_max)
    @constraint(m, fixed[i in 1:n_generators], α[i] == alpha[i])
    #@constraint(m, reserve, sum(α[i] for i in 1:n_generators) == 1)
    @constraint(m, cc1[i in 1:n_generators], p[i] + z * α[i] * s <= generators[i].g_max)
    @constraint(m, cc2[i in 1:n_generators], -p[i] + z * α[i] * s <= 0)

    # Linear Cost
    @expression(m, linear_cost, sum(p[i] * generators[i].cost for i in 1:n_generators))

    # Uncertain Cost
    @variable(m, r_uncert >= 0)
    @variable(m, r_sched >= 0)
    @constraint(m, vec(vcat(r_uncert, 0.5, C_rt' * α * s)) in RotatedSecondOrderCone())
    @constraint(m, vcat(r_sched, 0.5, C_rt' * p) in RotatedSecondOrderCone())
    @expression(m, quad_cost, r_sched + r_uncert)

    @objective(m, Min, linear_cost + quad_cost);

    return m

end
