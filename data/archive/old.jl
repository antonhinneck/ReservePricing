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
