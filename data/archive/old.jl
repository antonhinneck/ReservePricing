    #@constraint(m, vcat(0.5, r_sched, C_rt' * p_uncert) in RotatedSecondOrderCone())
    #@constraint(m, vec(vcat(r_uncert, C_rt * α * Σ_sq)) in RotatedSecondOrderCone())
    #@constraint(m, cost, vec(vcat(0.5, r_uncert, C_rt' * p_uncert)) in RotatedSecondOrderCone())


    @variable(m, p_uncert[1:n_generators] >= 0)
    #@constraint(m, conicgen[i in 1:n_generators], vec(vcat(p_uncert[i], α[i, :]' * Σ_sq)) in SecondOrderCone())
    #@constraint(m, feas, sum(p_uncert[i] for i in 1:n_generators) == s)

    #.+ pu_by_bus

    MSK_IPAR_INTPNT_SCALING = 1,
    MSK_IPAR_INTPNT_MAX_NUM_REFINEMENT_STEPS = 10, MSK_DPAR_INTPNT_TOL_STEP_SIZE = 0.999,
    MSK_DPAR_INTPNT_TOL_DSAFE = 2,
    #MSK_DPAR_INTPNT_TOL_DFEAS = 0.00001,
    #MSK_DPAR_INTPNT_CO_TOL_DFEAS = 0.00001,
    #MSK_DPAR_INTPNT_CO_TOL_INFEAS = 0.000072,
    #MSK_DPAR_INTPNT_CO_TOL_MU_RED = 1,
    #MSK_DPAR_INTPNT_CO_TOL_NEAR_REL = 10000,
    #MSK_DPAR_INTPNT_CO_TOL_PFEAS = 0.00001,
    #MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 0.00001,
    MSK_DPAR_INTPNT_TOL_REL_STEP = 0.9,
    MSK_IPAR_INTPNT_DIFF_STEP = 1

    #@constraint(m, uncert_gen[i in 1:n_generators], p_uncert[i] == sum(α[i, u] * Σ[u,u] for u in 1:n_farms))
    #@constraint(m, uncert_gen1[i in 1:n_generators], pp_uncert[i] == sum(αp[i, u] * Σ[u,u] for u in 1:n_farms))
    #@constraint(m, uncert_gen2[i in 1:n_generators], pm_uncert[i] == sum(αm[i, u] * Σ[u,u] for u in 1:n_farms))

    #=
    # Create unvertain generation (wind farms)
    wp = 1.25
    factor_σ =  1.25 * wp

    farms = []
    #                 capa    mva                   var    mva  bus
    #                         Base                         Base
    push!(farms, Farm(70.0  / 100 * wp, factor_σ * 7.0  / 100, 3))
    push!(farms, Farm(147.0 / 100 * wp, factor_σ * 14.7 / 100, 8))
    push!(farms, Farm(102.0 / 100 * wp, factor_σ * 10.2 / 100, 11))
    push!(farms, Farm(105.0 / 100 * wp, factor_σ * 10.5 / 100, 20))
    push!(farms, Farm(113.0 / 100 * wp, factor_σ * 11.3 / 100, 24))
    push!(farms, Farm(84.0  / 100 * wp, factor_σ * 8.4  / 100, 26))
    push!(farms, Farm(59.0  / 100 * wp, factor_σ * 5.9  / 100, 31))
    push!(farms, Farm(250.0 / 100 * wp, factor_σ * 25.0 / 100, 38))
    push!(farms, Farm(118.0 / 100 * wp, factor_σ * 11.8 / 100, 43))
    push!(farms, Farm(76.0  / 100 * wp, factor_σ * 7.6  / 100, 49))
    push!(farms, Farm(72.0  / 100 * wp, factor_σ * 7.2  / 100, 53))=#
