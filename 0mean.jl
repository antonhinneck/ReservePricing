include("models/0mean/dccc.jl")
m_dccc = build_dccc(generators, buses, lines, uRESs)
optimize!(m_dccc)
termination_status(m_dccc)
z1 = objective_value(m_dccc)
value(m_dccc[:det_c])
value(m_dccc[:d_lin])
value(m_dccc[:d_con])
value(m_dccc[:d_quad])
value(m_dccc[:unc_c])
χ0 = dual.(m_dccc[:γ])
δ1ms = dual.(m_dccc[:cc1])
δ1ps = dual.(m_dccc[:cc2])
sum(δ1ms)
sum(δ1ps)

include("models/0mean/dccc_apx.jl")
mylim = 0.0000
m_dccc_apx = build_dccc_apx(generators, buses, lines, uRESs, α_min = [mylim for i in 1:n_generators])
optimize!(m_dccc_apx)
termination_status(m_dccc_apx)
objective_value(m_dccc_apx)
sum(value.(m_dccc_apx[:p]))
value.(m_dccc_apx[:α])
sum(value.(m_dccc_apx[:α]))
sum(value.(m_dccc_apx[:ψ]))
value(m_dccc_apx[:det_c])
value(m_dccc_apx[:d_lin])
value(m_dccc_apx[:d_con])
value(m_dccc_apx[:d_quad])
value(m_dccc_apx[:unc_c])
value(m_dccc_apx[:d_bil])
χ1 = dual.(m_dccc_apx[:γ])

include("models/0mean/dccc_det.jl")
m_dccc_det = build_dccc_det(generators, buses, lines, uRESs, value.(m_dccc_apx[:p]))
optimize!(m_dccc_det)
termination_status(m_dccc_det)
objective_value(m_dccc_det)
value(m_dccc_det[:det_c])
value(m_dccc_det[:d_lin])
value(m_dccc_det[:d_con])
value(m_dccc_det[:d_quad])
value(m_dccc_det[:unc_c])
value(m_dccc_det[:d_bil])
χ1 = dual.(m_dccc_det[:γ])
δ1ms = dual.(m_dccc_det[:cc1])
δ1ps = dual.(m_dccc_det[:cc2])
sum(δ1ms)
sum(δ1ps)

include("models/0mean/dccc_a_apx.jl")
m_dccc_a_apx = build_dccc_a_apx(generators, buses, lines, uRESs, αm_min = zeros(n_generators), αp_min = zeros(n_generators))
optimize!(m_dccc_a_apx)
termination_status(m_dccc_a_apx)
z4 = objective_value(m_dccc_a_apx)
value(m_dccc_a_apx[:det_c])
value(m_dccc_a_apx[:d_lin])
value(m_dccc_a_apx[:d_con])
value(m_dccc_a_apx[:d_quad])
value(m_dccc_a_apx[:unc_c])
value(m_dccc_a_apx[:d_bil])
value.(m_dccc_a_apx[:αm])
value.(m_dccc_a_apx[:αp])
sum(value.(m_dccc_a_apx[:p]))
sum(value.(m_dccc_det[:p]))
value.(m_dccc_a_apx[:p])
d2 = dual.(m_dccc_a_apx[:γm])
d2 = dual.(m_dccc_a_apx[:γp])
value.(m_dccc_a_apx[:p])

include("models/0mean/dccc_a_det.jl")
m_dccc_a_det = build_dccc_a_det(generators, buses, lines, uRESs, value.(m_dccc_a_apx[:p]))
optimize!(m_dccc_a_det)
termination_status(m_dccc_a_det)
z4 = objective_value(m_dccc_a_det)
value(m_dccc_a_det[:det_c])
value(m_dccc_a_det[:d_lin])
value(m_dccc_a_det[:d_con])
value(m_dccc_a_det[:d_quad])
value(m_dccc_a_det[:unc_c])
value(m_dccc_a_det[:d_bil])
χm = dual.(m_dccc_a_det[:γm])
χp = dual.(m_dccc_a_det[:γp])

value.(m_dccc_det[:α])
value.(m_dccc_a_det[:αp])

function sw2n2n(alpha)

    # g x f
    new_alpha = ones(n_generators, length(uRESs))
    for i in 1:length(uRESs)
        new_alpha[:, i] = [alpha[i] for i in 1:n_generators]
    end
    return new_alpha

end

include("models/0mean/dccc_n2n.jl")
m_dccc_n2n = build_dccc_n2n(generators, buses, lines, uRESs)
optimize!(m_dccc_n2n)
termination_status(m_dccc_n2n)
z4 = objective_value(m_dccc_n2n)
value(m_dccc_n2n[:det_c])
value(m_dccc_n2n[:d_lin])
value(m_dccc_n2n[:d_con])
value(m_dccc_n2n[:d_quad])
value(m_dccc_n2n[:unc_c])
value.(m_dccc_n2n[:p_uncert])
χ0u = dual.(m_dccc_n2n[:χ])
δ0m = dual.(m_dccc_n2n[:cc1])
δ0p = dual.(m_dccc_n2n[:cc2])
sum(δm)
sum(δp)


sum(χ1u)
for i in 1:11
    print(string(round(χ1u[i], digits = 1)," & "))
end

α_min_init = zeros((n_generators, length(uRESs))) .+ 0.001 # sw2n2n(value.(m_dccc_a_det[:αp])) .- 0.000 #ones((n_generators, n_farms)) * 0.000001
α_max_init = ones((n_generators, length(uRESs)))
include("models/0mean/dccc_n2n_apx.jl")
m_dccc_n2n_apx = build_dccc_n2n_apx(generators, buses, lines, uRESs, α_min_init, α_max_init)
optimize!(m_dccc_n2n_apx)
termination_status(m_dccc_n2n_apx)
z4 = objective_value(m_dccc_n2n_apx)
value(m_dccc_n2n_apx[:det_c])
value(m_dccc_n2n_apx[:d_lin])
value(m_dccc_n2n_apx[:d_con])
value(m_dccc_n2n_apx[:d_quad])
value(m_dccc_n2n_apx[:unc_c])
value(m_dccc_n2n_apx[:u_bil])
value.(m_dccc_n2n_apx[:p_uncert])
d31 = dual.(m_dccc_n2n_apx[:χ])

α_det = value.(m_dccc_n2n_apx[:α])
p_det = value.(m_dccc_n2n_apx[:p])

include("models/0mean/dccc_n2n_det.jl")
m_dccc_n2n_det = build_dccc_n2n_det(generators, buses, lines, uRESs, α_det, p_det = p_det)
optimize!(m_dccc_n2n_det )
termination_status(m_dccc_n2n_det)
z4 = objective_value(m_dccc_n2n_det)
value(m_dccc_n2n_det[:det_c])
value(m_dccc_n2n_det[:d_lin])
value(m_dccc_n2n_det[:d_con])
value(m_dccc_n2n_det[:d_quad])
value(m_dccc_n2n_det[:unc_c])
value(m_dccc_n2n_det[:u_bil])
value(m_dccc_n2n_det[:u_quad])
χ1u = dual.(m_dccc_n2n_det[:χ])
δ1m = dual.(m_dccc_n2n_det[:cc1])
δ1p = dual.(m_dccc_n2n_det[:cc2])
sum(δ1m)
sum(δ1p)

α_min_initm = zeros((n_generators, n_farms)) .+ 0.0001 # value.(m_dccc_n2n_det[:α]) # ones((n_generators, n_farms)) * 0.01
α_max_initm = ones((n_generators, n_farms)) * 1.0
α_min_initp = zeros((n_generators, n_farms)) .+ 0.0001 # value.(m_dccc_n2n_det[:α]) # ones((n_generators, n_farms)) * 0.01
α_max_initp = ones((n_generators, n_farms)) * 1.0

include("models/0mean/dccc_n2n_a_apx_alpha.jl")
m_dccc_n2n_a_apx_alpha = build_dccc_n2n_a_apx_alpha(generators, buses, lines, uRESs, α_min_initm, α_max_initm, α_min_initp, α_max_initp)
optimize!(m_dccc_n2n_a_apx_alpha)
objective_value(m_dccc_n2n_a_apx_alpha)
termination_status(m_dccc_n2n_a_apx_alpha)
value(m_dccc_n2n_a_apx_alpha[:det_c])
value(m_dccc_n2n_a_apx_alpha[:d_lin])
value(m_dccc_n2n_a_apx_alpha[:d_con])
value(m_dccc_n2n_a_apx_alpha[:d_quad])
value(m_dccc_n2n_a_apx_alpha[:unc_c])
value(m_dccc_n2n_a_apx_alpha[:u_bil])
d41 = dual.(m_dccc_n2n_a_apx_alpha[:χp])
d42 = dual.(m_dccc_n2n_a_apx_alpha[:χm])

α_detm = value.(m_dccc_n2n_a_apx_alpha[:αm])
α_detp = value.(m_dccc_n2n_a_apx_alpha[:αp])
p_det = value.(m_dccc_n2n_a_apx_alpha[:p])
# p_det = value.(m_dccc_n2n_det[:p])

include("models/0mean/dccc_n2n_a_det_p.jl")
m_dccc_n2n_a_det_p = build_dccc_n2n_a_det_p(generators, buses, lines, uRESs, p_det)
optimize!(m_dccc_n2n_a_det_p)
termination_status(m_dccc_n2n_a_det_p)
objective_value(m_dccc_n2n_a_det_p)
value(m_dccc_n2n_a_det_p[:det_c])
value(m_dccc_n2n_a_det_p[:d_lin])
value(m_dccc_n2n_a_det_p[:d_con])
value(m_dccc_n2n_a_det_p[:d_quad])
value(m_dccc_n2n_a_det_p[:unc_c])
value(m_dccc_n2n_a_det_p[:u_quad])
value(m_dccc_n2n_a_det_p[:u_bil])
χ1up = dual.(m_dccc_n2n_a_det_p[:χp])
χ1um = dual.(m_dccc_n2n_a_det_p[:χm])
