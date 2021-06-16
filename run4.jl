ss = 0.5
include("load_data.jl")

include("models/dccc.jl")
m_dccc_cheb = build_dccc_cheb(generators, buses, lines, uRESs)
optimize!(m_dccc_cheb)
objective_value(m_dccc_cheb)
γc3 = dual.(m_dccc_cheb[:γ])
lmps13 = dual.(m_dccc_cheb[:mc])
cc1_13 = dual.(m_dccc_cheb[:cc1])
cc2_13 = dual.(m_dccc_cheb[:cc2])

include("models/dccc_a_apx.jl")
m_dccc_a_apx = build_dccc_a_apx(generators, buses, lines, uRESs, αm_min = zeros(n_generators), αp_min = zeros(n_generators), cheb = true)
optimize!(m_dccc_a_apx)
objective_value(m_dccc_a_apx)
