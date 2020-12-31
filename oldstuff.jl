# 
# function adjust_limits_a(alphas::Array{T, 2} where T <: Real; lim = Float64(0.0), dec = 6)
#
#     new_min = deepcopy(alphas)
#     new_max = deepcopy(alphas)
#
#     new_min -= ones(size(alphas)) * lim
#     new_max += ones(size(alphas)) * lim
#
#     return new_min, new_max
# end
#
# function min_alg_a(iters::Int64, α_initm, α_initp; lim_init = 0.01, lim_const = 0.01, scheduler = 1.0, verbose = true)
#
#     progression_apx = Vector{Float64}()
#     progression_det = Vector{Float64}()
#     iteration = Vector{String}()
#
#     # if α_initm == nothing
#     #     α_initm = zeros(n_generators, n_farms)
#     #     α_initp = zeros(n_generators, n_farms)
#     # end
#
#
#     lims_min_m = nothing
#     lims_max_m = nothing
#
#     lims_min_p = nothing
#     lims_max_p = nothing
#
#     current_apx_model = nothing
#     ctr = 0
#
#     for i in 1:iters
#
#         if i == 1
#             #m_dccc_n2n_a_apx_alg = build_dccc_n2n_a_apx(generators, buses, lines, farms, ones(n_generators, n_farms) * lim_init, ones(n_generators, n_farms) * 1.0, ones(n_generators, n_farms) * lim_init, ones(n_generators, n_farms) * 1.0
#             m_dccc_n2n_a_apx_alg = build_dccc_n2n_a_apx(generators, buses, lines, farms, α_initm .- lim_init, ones(n_generators, n_farms), α_initp .- lim_init, ones(n_generators, n_farms))
#             optimize!(m_dccc_n2n_a_apx_alg)
#             println(objective_value(m_dccc_n2n_a_apx_alg))
#         else
#             m_dccc_n2n_a_apx_alg = build_dccc_n2n_a_apx(generators, buses, lines, farms, lims_min_m, lims_max_m, lims_min_p, lims_max_p)
#             optimize!(m_dccc_n2n_a_apx_alg)
#         end
#
#         m_dccc_n2n_a_det_alg = build_dccc_n2n_a_det(generators, buses, lines, farms, value.(m_dccc_n2n_a_apx_alg[:αm]), value.(m_dccc_n2n_a_apx_alg[:αp]), output_level = 0)
#         optimize!(m_dccc_n2n_a_det_alg)
#
#         if verbose
#             println(string("ITR: ", i, " APX: ", objective_value(m_dccc_n2n_a_apx_alg), " DET: ", objective_value(m_dccc_n2n_a_det_alg)))
#         end
#
#         lims_min_m, lims_max_m = adjust_limits_a(value.(m_dccc_n2n_a_apx_alg[:αm]), lim = lim_const)
#         lims_min_p, lims_max_p = adjust_limits_a(value.(m_dccc_n2n_a_apx_alg[:αp]), lim = lim_const)
#
#         if i == 1
#             push!(progression_apx, objective_value(m_dccc_n2n_a_apx_alg))
#             push!(progression_det, objective_value(m_dccc_n2n_a_apx_alg))
#             push!(iteration, string(1))
#             ctr += 1
#         else
#             if progression_det[ctr] > objective_value(m_dccc_n2n_a_det_alg) && termination_status(m_dccc_n2n_a_det_alg) == MOI.TerminationStatusCode(1) && objective_value(m_dccc_n2n_a_det_alg) != 0.0
#                 push!(iteration, string(i))
#                 push!(progression_apx, objective_value(m_dccc_n2n_a_apx_alg))
#                 push!(progression_det, objective_value(m_dccc_n2n_a_apx_alg))
#                 current_apx_model = m_dccc_n2n_a_apx_alg
#                 ctr += 1
#             end
#         end
#
#         #lim_const = lim_const * scheduler
#     end
#
#     return current_apx_model, [progression_apx, progression_det, iteration]
# end

# value.(m_dccc_det[:α])
# value.(m_dccc_a_det[:αp])
#
# function sw2n2n(alpha)
#
#     # g x f
#     new_alpha = ones(n_generators, n_farms)
#     for i in 1:n_farms
#         new_alpha[:, i] = [alpha[i] for i in 1:n_generators]
#     end
#     return new_alpha
#
# end
#
# function adjust_limits_sw(alphas::Array{T, 1} where T <: Real; lim = Float64(0.0), dec = 6)
#
#     new_min = zeros(length(alphas))
#     new_max = zeros(length(alphas))
#
#     for (i, a) in enumerate(alphas)
#         new_min[i] = max(round(a - lim, digits = dec), 0.0)
#         new_max[i] = min(round(a + lim, digits = dec), 1.0)
#     end
#
#     return new_min, new_max
# end
#
# function min_alg_sws(iters::Int64; lim_init = 0.00001, lim_const = 0.001, scheduler = 1.0, verbose = true)
#
#     progression_apx = Vector{Float64}()
#     progression_det = Vector{Float64}()
#     iteration = Vector{String}()
#
#     lims_min = nothing
#     lims_max = nothing
#
#     current_apx_model = nothing
#     ctr = 0
#
#     for i in 1:iters
#
#         if i == 1
#             m_dccc_apx_alg = build_dccc_apx(generators, buses, lines, farms, α_min = [lim_init for i in 1:n_generators])
#             optimize!(m_dccc_apx_alg)
#         else
#             m_dccc_apx_alg = build_dccc_apx(generators, buses, lines, farms, α_min = lims_min, α_max = lims_max, output_level = 0)
#             optimize!(m_dccc_apx_alg)
#         end
#
#         m_dccc_det_alg = build_dccc_det(generators, buses, lines, farms, value.(m_dccc_apx_alg[:α]), value.(m_dccc_apx_alg[:p]), output_level = 0)
#         optimize!(m_dccc_det_alg)
#
#         if verbose
#             println(string("ITR: ", i, " APX: ", objective_value(m_dccc_apx_alg), " DET: ", objective_value(m_dccc_det_alg)))
#         end
#
#         lims_min, lims_max = adjust_limits_sw(value.(m_dccc_apx_alg[:α]), lim = lim_const)
#
#         if i == 1
#             push!(progression_apx, objective_value(m_dccc_apx_alg))
#             push!(progression_det, objective_value(m_dccc_det_alg))
#             push!(iteration, string(1))
#             ctr += 1
#         else
#             if progression_det[ctr] > objective_value(m_dccc_det_alg) && termination_status(m_dccc_det_alg) == MOI.TerminationStatusCode(1) && objective_value(m_dccc_det_alg) != 0.0
#                 push!(iteration, string(i))
#                 push!(progression_apx, objective_value(m_dccc_apx_alg))
#                 push!(progression_det, objective_value(m_dccc_det_alg))
#                 current_apx_model = m_dccc_apx_alg
#                 ctr += 1
#             end
#         end
#
#         #lim_const = lim_const * scheduler
#     end
#
#     return current_apx_model, [progression_apx, progression_det, iteration]
# end

# min_model, data = min_alg_sws(60, lim_init = 0.00, lim_const = 0.006)
# objective_value(min_model)

# m_dccc_det = build_dccc_det(generators, buses, lines, farms, value.(min_model[:α]), value.(min_model[:p]), output_level = 1)
# optimize!(m_dccc_det)
# objective_value(m_dccc_det)
# value(m_dccc_det[:det_c])
# value(m_dccc_det[:d_lin])
# value(m_dccc_det[:d_con])
# value(m_dccc_det[:d_quad])
# value(m_dccc_det[:unc_c])
# value(m_dccc_det[:d_bil])

using PyPlot

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.09,0.2,0.905,0.78])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("Iteration")
ylabel("\$z\$")
ylim(bottom = 98260, top = 98550)
#ylabel("\$\\chi^{+}_{u}\$")

plot(data[3], data[1], color = "lightgreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$z^{APX}\$")
plot(data[3], data[2], color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$z^{DET}\$")

legend(loc = "upper right", fancybox = false, edgecolor = "black", framealpha = 0.9)
savefig(string("plots_final//algorithm.pdf"), format = :pdf)
#sm = sqrt(ζ2) * s
