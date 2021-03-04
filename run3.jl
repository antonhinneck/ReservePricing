ss = 0.5
include("load_data.jl")

α_min_initm = zeros((n_generators, n_farms)) .+ 0.0
α_max_initm = ones((n_generators, n_farms)) * 1.0
α_min_initp = zeros((n_generators, n_farms)) .+ 0.0
α_max_initp = ones((n_generators, n_farms)) * 1.0

##
##
##

include("models/dccc_cheb.jl")
m_dccc_cheb = build_dccc_cheb(generators, buses, lines, uRESs)
optimize!(m_dccc_cheb)
γc1 = dual.(m_dccc_cheb[:γ])
lmps1 = dual.(m_dccc_cheb[:mc])
cc1_11 = dual.(m_dccc_cheb[:cc1])
cc2_11 = dual.(m_dccc_cheb[:cc2])

include("models/dccc_n2n_cheb.jl")
m_dccc_n2n_cheb = build_dccc_n2n_cheb(generators, buses, lines, uRESs)
optimize!(m_dccc_n2n_cheb)
lmps3 = dual.(m_dccc_cheb[:mc])
χc1 = dual.(m_dccc_n2n_cheb[:χ])
cc1_31 = dual.(m_dccc_n2n_cheb[:cc1])
cc2_31 = dual.(m_dccc_n2n_cheb[:cc2])

include("models/dccc_a_apx.jl")
m_dccc_a_apx = build_dccc_a_apx(generators, buses, lines, uRESs, αm_min = zeros(n_generators), αp_min = zeros(n_generators), cheb = true)
optimize!(m_dccc_a_apx)

include("models/dccc_a_det_cheb.jl")
m_dccc_a_det_cheb = build_dccc_a_det_cheb(generators, buses, lines, uRESs, value.(m_dccc_a_apx[:p]))
optimize!(m_dccc_a_det_cheb)
γmc1 = dual.(m_dccc_a_det_cheb[:γm])
γpc1 = dual.(m_dccc_a_det_cheb[:γp])

include("models/dccc_a_det_falpha.jl")
m_dccc_a_det_falpha = build_dccc_a_det_falpha(generators, buses, lines, uRESs; apdet = value.(m_dccc_a_det_cheb[:αp]), amdet = value.(m_dccc_a_det_cheb[:αm]), cheb = true)
optimize!(m_dccc_a_det_falpha)
z2c = objective_value(m_dccc_a_det_falpha)
lmps21 = dual.(m_dccc_a_det_falpha[:mc])
cc1_21 = dual.(m_dccc_n2n_cheb[:cc1])
cc2_21 = dual.(m_dccc_n2n_cheb[:cc2])

include("models/dccc_n2n_a_apx_alpha.jl")
m_dccc_n2n_a_apx_alpha = build_dccc_n2n_a_apx_alpha(generators, buses, lines, uRESs, α_min_initm, α_max_initm, α_min_initp, α_max_initp, cheb = true)
optimize!(m_dccc_n2n_a_apx_alpha)

p_det = value.(m_dccc_n2n_a_apx_alpha[:p])

include("models/dccc_n2n_a_det_p_cheb.jl")
m_dccc_n2n_a_det_p_cheb = build_dccc_n2n_a_det_p_cheb(generators, buses, lines, uRESs, value.(m_dccc_n2n_a_apx_alpha[:p]))
optimize!(m_dccc_n2n_a_det_p_cheb)
χ1up_cheb1 = dual.(m_dccc_n2n_a_det_p_cheb[:χp])
χ1um_cheb1 = dual.(m_dccc_n2n_a_det_p_cheb[:χm])

include("models/dccc_n2n_a_det_p_falpha.jl")
m_dccc_n2n_a_det_p_falpha = build_dccc_n2n_a_det_p_fa(generators, buses, lines, uRESs, abs.(value.(m_dccc_n2n_a_det_p_cheb[:αm])), abs.(value.(m_dccc_n2n_a_det_p_cheb[:αp])), cheb = true)
optimize!(m_dccc_n2n_a_det_p_falpha)
lmps41 = dual.(m_dccc_n2n_a_det_p_falpha[:mc])
cc1_41 = dual.(m_dccc_n2n_a_det_p_falpha[:cc1])
cc2_41 = dual.(m_dccc_n2n_a_det_p_falpha[:cc2])

##
##
##

ss = 1.0
include("load_data.jl")

include("models/dccc_cheb.jl")
m_dccc_cheb = build_dccc_cheb(generators, buses, lines, uRESs)
optimize!(m_dccc_cheb)
γc2 = dual.(m_dccc_cheb[:γ])
lmps12 = dual.(m_dccc_cheb[:mc])
cc1_12 = dual.(m_dccc_cheb[:cc1])
cc2_12 = dual.(m_dccc_cheb[:cc2])

include("models/dccc_n2n_cheb.jl")
m_dccc_n2n_cheb = build_dccc_n2n_cheb(generators, buses, lines, uRESs)
optimize!(m_dccc_n2n_cheb)
lmps32 = dual.(m_dccc_cheb[:mc])
χc2 = dual.(m_dccc_n2n_cheb[:χ])
cc1_32 = dual.(m_dccc_n2n_cheb[:cc1])
cc2_32 = dual.(m_dccc_n2n_cheb[:cc2])

include("models/dccc_a_apx.jl")
m_dccc_a_apx = build_dccc_a_apx(generators, buses, lines, uRESs, αm_min = zeros(n_generators), αp_min = zeros(n_generators), cheb = true)
optimize!(m_dccc_a_apx)

include("models/dccc_a_det_cheb.jl")
m_dccc_a_det_cheb = build_dccc_a_det_cheb(generators, buses, lines, uRESs, value.(m_dccc_a_apx[:p]))
optimize!(m_dccc_a_det_cheb)
γmc2 = dual.(m_dccc_a_det_cheb[:γm])
γpc2 = dual.(m_dccc_a_det_cheb[:γp])

include("models/dccc_a_det_falpha.jl")
m_dccc_a_det_falpha = build_dccc_a_det_falpha(generators, buses, lines, uRESs; apdet = value.(m_dccc_a_det_cheb[:αp]), amdet = value.(m_dccc_a_det_cheb[:αm]), cheb = true)
optimize!(m_dccc_a_det_falpha)
z2c = objective_value(m_dccc_a_det_falpha)
lmps22 = dual.(m_dccc_a_det_falpha[:mc])
cc1_22 = dual.(m_dccc_n2n_cheb[:cc1])
cc2_22 = dual.(m_dccc_n2n_cheb[:cc2])

include("models/dccc_n2n_a_apx_alpha.jl")
m_dccc_n2n_a_apx_alpha = build_dccc_n2n_a_apx_alpha(generators, buses, lines, uRESs, α_min_initm, α_max_initm, α_min_initp, α_max_initp, cheb = true)
optimize!(m_dccc_n2n_a_apx_alpha)

p_det = value.(m_dccc_n2n_a_apx_alpha[:p])

include("models/dccc_n2n_a_det_p_cheb.jl")
m_dccc_n2n_a_det_p_cheb = build_dccc_n2n_a_det_p_cheb(generators, buses, lines, uRESs, value.(m_dccc_n2n_a_apx_alpha[:p]))
optimize!(m_dccc_n2n_a_det_p_cheb)
χ1up_cheb1 = dual.(m_dccc_n2n_a_det_p_cheb[:χp])
χ1um_cheb1 = dual.(m_dccc_n2n_a_det_p_cheb[:χm])

include("models/dccc_n2n_a_det_p_falpha.jl")
m_dccc_n2n_a_det_p_falpha = build_dccc_n2n_a_det_p_fa(generators, buses, lines, uRESs, abs.(value.(m_dccc_n2n_a_det_p_cheb[:αm])), abs.(value.(m_dccc_n2n_a_det_p_cheb[:αp])), cheb = true)
optimize!(m_dccc_n2n_a_det_p_falpha)
lmps42 = dual.(m_dccc_n2n_a_det_p_falpha[:mc])
cc1_42 = dual.(m_dccc_n2n_a_det_p_falpha[:cc1])
cc2_42 = dual.(m_dccc_n2n_a_det_p_falpha[:cc2])

##
##
##

ss = 2.0
include("load_data.jl")

include("models/dccc_cheb.jl")
m_dccc_cheb = build_dccc_cheb(generators, buses, lines, uRESs)
optimize!(m_dccc_cheb)
γc3 = dual.(m_dccc_cheb[:γ])
lmps13 = dual.(m_dccc_cheb[:mc])
cc1_13 = dual.(m_dccc_cheb[:cc1])
cc2_13 = dual.(m_dccc_cheb[:cc2])

include("models/dccc_n2n_cheb.jl")
m_dccc_n2n_cheb = build_dccc_n2n_cheb(generators, buses, lines, uRESs)
optimize!(m_dccc_n2n_cheb)
lmps33 = dual.(m_dccc_cheb[:mc])
χc3 = dual.(m_dccc_n2n_cheb[:χ])
cc1_33 = dual.(m_dccc_n2n_cheb[:cc1])
cc2_33 = dual.(m_dccc_n2n_cheb[:cc2])

include("models/dccc_a_apx.jl")
m_dccc_a_apx = build_dccc_a_apx(generators, buses, lines, uRESs, αm_min = zeros(n_generators), αp_min = zeros(n_generators), cheb = true)
optimize!(m_dccc_a_apx)

include("models/dccc_a_det_cheb.jl")
m_dccc_a_det_cheb = build_dccc_a_det_cheb(generators, buses, lines, uRESs, value.(m_dccc_a_apx[:p]))
optimize!(m_dccc_a_det_cheb)
γmc3 = dual.(m_dccc_a_det_cheb[:γm])
γpc3 = dual.(m_dccc_a_det_cheb[:γp])

include("models/dccc_a_det_falpha.jl")
m_dccc_a_det_falpha = build_dccc_a_det_falpha(generators, buses, lines, uRESs; apdet = value.(m_dccc_a_det_cheb[:αp]), amdet = value.(m_dccc_a_det_cheb[:αm]), cheb = true)
optimize!(m_dccc_a_det_falpha)
z2c = objective_value(m_dccc_a_det_falpha)
lmps23 = dual.(m_dccc_a_det_falpha[:mc])
cc1_23 = dual.(m_dccc_n2n_cheb[:cc1])
cc2_23 = dual.(m_dccc_n2n_cheb[:cc2])

include("models/dccc_n2n_a_apx_alpha.jl")
m_dccc_n2n_a_apx_alpha = build_dccc_n2n_a_apx_alpha(generators, buses, lines, uRESs, α_min_initm, α_max_initm, α_min_initp, α_max_initp, cheb = true)
optimize!(m_dccc_n2n_a_apx_alpha)

p_det = value.(m_dccc_n2n_a_apx_alpha[:p])

include("models/dccc_n2n_a_det_p_cheb.jl")
m_dccc_n2n_a_det_p_cheb = build_dccc_n2n_a_det_p_cheb(generators, buses, lines, uRESs, value.(m_dccc_n2n_a_apx_alpha[:p]))
optimize!(m_dccc_n2n_a_det_p_cheb)
χ1up_cheb3 = dual.(m_dccc_n2n_a_det_p_cheb[:χp])
χ1um_cheb3 = dual.(m_dccc_n2n_a_det_p_cheb[:χm])

include("models/dccc_n2n_a_det_p_falpha.jl")
m_dccc_n2n_a_det_p_falpha = build_dccc_n2n_a_det_p_fa(generators, buses, lines, uRESs, abs.(value.(m_dccc_n2n_a_det_p_cheb[:αm])), abs.(value.(m_dccc_n2n_a_det_p_cheb[:αp])), cheb = true)
optimize!(m_dccc_n2n_a_det_p_falpha)
lmps43 = dual.(m_dccc_n2n_a_det_p_falpha[:mc])
cc1_43 = dual.(m_dccc_n2n_a_det_p_falpha[:cc1])
cc2_43 = dual.(m_dccc_n2n_a_det_p_falpha[:cc2])

##
##
##

ss = 4.0
include("load_data.jl")

include("models/dccc_cheb.jl")
m_dccc_cheb = build_dccc_cheb(generators, buses, lines, uRESs)
optimize!(m_dccc_cheb)
γc4 = dual.(m_dccc_cheb[:γ])
lmps14 = dual.(m_dccc_cheb[:mc])
cc1_14 = dual.(m_dccc_cheb[:cc1])
cc2_14 = dual.(m_dccc_cheb[:cc2])

include("models/dccc_n2n_cheb.jl")
m_dccc_n2n_cheb = build_dccc_n2n_cheb(generators, buses, lines, uRESs)
optimize!(m_dccc_n2n_cheb)
lmps34 = dual.(m_dccc_cheb[:mc])
χc4 = dual.(m_dccc_n2n_cheb[:χ])
cc1_34 = dual.(m_dccc_n2n_cheb[:cc1])
cc2_34 = dual.(m_dccc_n2n_cheb[:cc2])

include("models/dccc_a_apx.jl")
m_dccc_a_apx = build_dccc_a_apx(generators, buses, lines, uRESs, αm_min = zeros(n_generators), αp_min = zeros(n_generators), cheb = true)
optimize!(m_dccc_a_apx)

include("models/dccc_a_det_cheb.jl")
m_dccc_a_det_cheb = build_dccc_a_det_cheb(generators, buses, lines, uRESs, value.(m_dccc_a_apx[:p]))
optimize!(m_dccc_a_det_cheb)
γmc4 = dual.(m_dccc_a_det_cheb[:γm])
γpc4 = dual.(m_dccc_a_det_cheb[:γp])

include("models/dccc_a_det_falpha.jl")
m_dccc_a_det_falpha = build_dccc_a_det_falpha(generators, buses, lines, uRESs; apdet = value.(m_dccc_a_det_cheb[:αp]), amdet = value.(m_dccc_a_det_cheb[:αm]), cheb = true)
optimize!(m_dccc_a_det_falpha)
z2c = objective_value(m_dccc_a_det_falpha)
lmps24 = dual.(m_dccc_a_det_falpha[:mc])
cc1_24 = dual.(m_dccc_n2n_cheb[:cc1])
cc2_24 = dual.(m_dccc_n2n_cheb[:cc2])

include("models/dccc_n2n_a_apx_alpha.jl")
m_dccc_n2n_a_apx_alpha = build_dccc_n2n_a_apx_alpha(generators, buses, lines, uRESs, α_min_initm, α_max_initm, α_min_initp, α_max_initp, cheb = true)
optimize!(m_dccc_n2n_a_apx_alpha)

p_det = value.(m_dccc_n2n_a_apx_alpha[:p])

include("models/dccc_n2n_a_det_p_cheb.jl")
m_dccc_n2n_a_det_p_cheb = build_dccc_n2n_a_det_p_cheb(generators, buses, lines, uRESs, value.(m_dccc_n2n_a_apx_alpha[:p]))
optimize!(m_dccc_n2n_a_det_p_cheb)
χ1up_cheb4 = dual.(m_dccc_n2n_a_det_p_cheb[:χp])
χ1um_cheb4 = dual.(m_dccc_n2n_a_det_p_cheb[:χm])

include("models/dccc_n2n_a_det_p_falpha.jl")
m_dccc_n2n_a_det_p_falpha = build_dccc_n2n_a_det_p_fa(generators, buses, lines, uRESs, abs.(value.(m_dccc_n2n_a_det_p_cheb[:αm])), abs.(value.(m_dccc_n2n_a_det_p_cheb[:αp])), cheb = true)
optimize!(m_dccc_n2n_a_det_p_falpha)
lmps44 = dual.(m_dccc_n2n_a_det_p_falpha[:mc])
cc1_44 = dual.(m_dccc_n2n_a_det_p_falpha[:cc1])
cc2_44 = dual.(m_dccc_n2n_a_det_p_falpha[:cc2])



fig = figure(figsize=(8, 2.8))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.085,0.15,0.91,0.82])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#yticks([0.001,0.01,0.1,0.0,1.0,1000.0], labels = ["\&10^{-3}\&","\&10^{-2}\&","\&10^{-1}\&","\&0\&","\&10^{1}\&","\&10^{3}\&"])
yticks([0.01, 0.0], labels = ["10", "0"])

ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$u\$")
ylim(bottom = 0.001, top = 35000000)
#xlim(left = 0.0, right = length(u_buses))
ylabel("\$[\\\$]\$")

plot(string(u_buses), χc1, color = "limegreen", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~38.8\\% \$")
plot(string(u_buses), [γc1 for i in u_buses], color = "limegreen", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~38.8\\% \$")

plot(string(u_buses), χc2, color = "silver", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~19.38\\% \$")
plot(string(u_buses), [γc2 for i in u_buses], color = "silver", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~19.38\\% \$")

plot(string(u_buses), χc3, color = "navy", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~9.7\\% \$")
plot(string(u_buses), [γc3 for i in u_buses], color = "navy", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~9.7\\% \$")

plot(string(u_buses), χc4, color = "crimson", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~4.85\\% \$")
plot(string(u_buses), [γc4 for i in u_buses], color = "crimson", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~4.85\\% \$")

# plot(string(u_buses), cc1, color = "limegreen", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~38.8\\% \$")
# plot(string(u_buses), [γc1 for i in u_buses], color = "limegreen", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~38.8\\% \$")

# plot(string(u_buses), χ1up, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi^{+}_{u}, Ncvx.\$")
# plot(string(u_buses), χ1um, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi^{-}_{u}, Ncvx.\$")

# plot(string(u_buses), χ1up_cheb, color = "lightgreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$Chebychev, up$")
# plot(string(u_buses), χ1um_cheb, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$Chebychev, down\$")
# plot(string(u_buses), χ1up, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\chi_{u}\$")
# plot(string(u_buses), χ1um, color = "cornflowerblue", lw = 1.2, ls = "dashed", label = "\$\\chi\$")

legend(loc = "upper center", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 4)
savefig(string("plots_final//prices_sym_scenarios1.pdf"), format = :pdf)

using PyPlot

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.09,0.22,0.905,0.74])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)
xticks([10,20,30,40,50], labels = ["10", "20", "40", "50"])

#yticks([0.001,0.01,0.1,0.0,1.0,1000.0], labels = ["\&10^{-3}\&","\&10^{-2}\&","\&10^{-1}\&","\&0\&","\&10^{1}\&","\&10^{3}\&"])
#yticks([0.01, 0.0], labels = ["10", "0"])

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$g\$")
ylim(bottom = 0.1, top = 1250)
#xlim(left = 0.0, right = length(u_buses))
ylabel("\$\\underline{\\delta}[\\\$]\$")

plot(ges, -cc2_11, color = "crimson", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~4.85\\% \$")
plot(ges, -cc2_12, color = "yellow", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~9.8\\% \$")
plot(ges, -cc2_13, color = "orange", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~19.76\\% \$")
plot(ges, -cc2_14, color = "silver", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~38.98\\% \$")

legend(loc = "upper center", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 4)
savefig(string("plots_final//delta_down_sym.pdf"), format = :pdf)

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.09,0.22,0.905,0.74])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)
xticks([10,20,30,40,50], labels = ["10", "20", "40", "50"])

#yticks([0.001,0.01,0.1,0.0,1.0,1000.0], labels = ["\&10^{-3}\&","\&10^{-2}\&","\&10^{-1}\&","\&0\&","\&10^{1}\&","\&10^{3}\&"])
#yticks([0.01, 0.0], labels = ["10", "0"])

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$g\$")
ylim(bottom = 0.1, top = 1250)
#xlim(left = 0.0, right = length(u_buses))
ylabel("\$\\bar{\\delta}[\\\$]\$")

plot(ges, -cc1_11, color = "crimson", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~4.85\\% \$")
plot(ges, -cc1_12, color = "yellow", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~9.8\\% \$")
plot(ges, -cc1_13, color = "orange", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~19.76\\% \$")
plot(ges, -cc1_14, color = "silver", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~38.98\\% \$")

legend(loc = "upper center", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 4)
savefig(string("plots_final//delta_up_sym.pdf"), format = :pdf)

##
##

fig = figure(figsize=(8, 2.8))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.085,0.15,0.91,0.82])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#yticks([0.001,0.01,0.1,0.0,1.0,1000.0], labels = ["\&10^{-3}\&","\&10^{-2}\&","\&10^{-1}\&","\&0\&","\&10^{1}\&","\&10^{3}\&"])
yticks([0.01, 0.0], labels = ["10", "0"])

ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$u\$")
ylim(bottom = 0.001, top = 60000000)
#xlim(left = 0.0, right = length(u_buses))
ylabel("\$[\\\$]\$")

plot(string(u_buses), -χ1um_cheb4, color = "crimson", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~4.85\\% \$")
plot(string(u_buses), [γmc4 for i in u_buses], color = "crimson", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~4.85\\% \$")

plot(string(u_buses), -χ1um_cheb3, color = "navy", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~9.7\\% \$")
plot(string(u_buses), [γmc3 for i in u_buses], color = "navy", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~9.7\\% \$")

plot(string(u_buses), -χ1um_cheb2, color = "silver", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~19.38\\% \$")
plot(string(u_buses), [γmc2 for i in u_buses], color = "silver", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~19.38\\% \$")

plot(string(u_buses), -χ1um_cheb1, color = "limegreen", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~38.8\\% \$")
plot(string(u_buses), [γmc1 for i in u_buses], color = "limegreen", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~38.8\\% \$")

# plot(string(u_buses), χ1up, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi^{+}_{u}, Ncvx.\$")
# plot(string(u_buses), χ1um, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi^{-}_{u}, Ncvx.\$")

# plot(string(u_buses), χ1up_cheb, color = "lightgreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$Chebychev, up$")
# plot(string(u_buses), χ1um_cheb, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$Chebychev, down\$")
# plot(string(u_buses), χ1up, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\chi_{u}\$")
# plot(string(u_buses), χ1um, color = "cornflowerblue", lw = 1.2, ls = "dashed", label = "\$\\chi\$")

legend(loc = "upper center", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 4)
savefig(string("plots_final//prices_sym_scenarios_minus.pdf"), format = :pdf)

##
##

fig = figure(figsize=(8, 2.8))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.085,0.15,0.91,0.82])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#yticks([0.001,0.01,0.1,0.0,1.0,1000.0], labels = ["\&10^{-3}\&","\&10^{-2}\&","\&10^{-1}\&","\&0\&","\&10^{1}\&","\&10^{3}\&"])
yticks([0.01, 0.0], labels = ["10", "0"])

ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$u\$")
ylim(bottom = 0.001, top = 60000000)
#xlim(left = 0.0, right = length(u_buses))
ylabel("\$[\\\$]\$")

plot(string(u_buses), χ1up_cheb4, color = "crimson", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~4.85\\% \$")
plot(string(u_buses), [γpc4 for i in u_buses], color = "crimson", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~4.85\\% \$")

plot(string(u_buses), χ1up_cheb3, color = "navy", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~9.7\\% \$")
plot(string(u_buses), [γpc3 for i in u_buses], color = "navy", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~9.7\\% \$")

plot(string(u_buses), χ1up_cheb2, color = "silver", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~19.38\\% \$")
plot(string(u_buses), [γpc2 for i in u_buses], color = "silver", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~19.38\\% \$")

plot(string(u_buses), χ1up_cheb1, color = "limegreen", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~38.8\\% \$")
plot(string(u_buses), [γpc1 for i in u_buses], color = "limegreen", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~38.8\\% \$")

# plot(string(u_buses), χ1up, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi^{+}_{u}, Ncvx.\$")
# plot(string(u_buses), χ1um, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi^{-}_{u}, Ncvx.\$")

# plot(string(u_buses), χ1up_cheb, color = "lightgreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$Chebychev, up$")
# plot(string(u_buses), χ1um_cheb, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$Chebychev, down\$")
# plot(string(u_buses), χ1up, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\chi_{u}\$")
# plot(string(u_buses), χ1um, color = "cornflowerblue", lw = 1.2, ls = "dashed", label = "\$\\chi\$")

legend(loc = "upper center", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 4)
savefig(string("plots_final//prices_sym_scenarios_plus.pdf"), format = :pdf)

##
##

fig = figure(figsize=(8, 2.8))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.085,0.15,0.91,0.82])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#yticks([0.001,0.01,0.1,0.0,1.0,1000.0], labels = ["\&10^{-3}\&","\&10^{-2}\&","\&10^{-1}\&","\&0\&","\&10^{1}\&","\&10^{3}\&"])
yticks([0.01, 0.0], labels = ["10", "0"])

ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$u\$")
ylim(bottom = 0.1, top = 60000000)
#xlim(left = 0.0, right = length(u_buses))
ylabel("\$[\\\$]\$")

plot(string(u_buses), χ1up_cheb4, color = "crimson", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~4.85\\% \$")
plot(string(u_buses), [γpc4 for i in u_buses], color = "crimson", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~4.85\\% \$")

plot(string(u_buses), χ1up_cheb3, color = "navy", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~9.7\\% \$")
plot(string(u_buses), [γpc3 for i in u_buses], color = "navy", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~9.7\\% \$")

plot(string(u_buses), χ1up_cheb2, color = "silver", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~19.38\\% \$")
plot(string(u_buses), [γpc2 for i in u_buses], color = "silver", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~19.38\\% \$")

plot(string(u_buses), χ1up_cheb1, color = "limegreen", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~38.8\\% \$")
plot(string(u_buses), [γpc1 for i in u_buses], color = "limegreen", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~38.8\\% \$")

# plot(string(u_buses), χ1up, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi^{+}_{u}, Ncvx.\$")
# plot(string(u_buses), χ1um, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi^{-}_{u}, Ncvx.\$")

# plot(string(u_buses), χ1up_cheb, color = "lightgreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$Chebychev, up$")
# plot(string(u_buses), χ1um_cheb, color = "lightseagreen", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$Chebychev, down\$")
# plot(string(u_buses), χ1up, color = "teal", lw = 1.2, ls = "dotted", marker = "D", ms = 4.0, mfc = "white", label = "\$\\chi_{u}\$")
# plot(string(u_buses), χ1um, color = "cornflowerblue", lw = 1.2, ls = "dashed", label = "\$\\chi\$")

legend(loc = "upper center", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 4)
savefig(string("plots_final//prices_sym_scenarios_plus.pdf"), format = :pdf)


##
##

using PyPlot
ges = [string(i) for i in 1:n_generators]

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.08,0.22,0.915,0.74])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)
xticks([10,20,30,40,50], labels = ["10", "20", "40", "50"])

#yticks([0.001,0.01,0.1,0.0,1.0,1000.0], labels = ["\&10^{-3}\&","\&10^{-2}\&","\&10^{-1}\&","\&0\&","\&10^{1}\&","\&10^{3}\&"])
#yticks([0.01, 0.0], labels = ["10", "0"])

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$g\$")
ylim(bottom = 0.1, top = 950)
#xlim(left = 0.0, right = length(u_buses))
ylabel("\$\\bar{\\delta}[\\\$]\$")

plot(ges, -cc1_41, color = "crimson", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~4.85\\% \$")
plot(ges, -cc1_42, color = "yellow", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~9.8\\% \$")
plot(ges, -cc1_43, color = "orange", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~19.76\\% \$")
plot(ges, -cc1_44, color = "silver", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~38.98\\% \$")
# plot(string(u_buses), [γpc4 for i in u_buses], color = "crimson", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~4.85\\% \$")
#
# plot(string(u_buses), χ1up_cheb3, color = "navy", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~9.7\\% \$")
# plot(string(u_buses), [γpc3 for i in u_buses], color = "navy", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~9.7\\% \$")
#
# plot(string(u_buses), χ1up_cheb2, color = "silver", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~19.38\\% \$")
# plot(string(u_buses), [γpc2 for i in u_buses], color = "silver", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~19.38\\% \$")
#
# plot(string(u_buses), χ1up_cheb1, color = "limegreen", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~38.8\\% \$")
# plot(string(u_buses), [γpc1 for i in u_buses], color = "limegreen", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~38.8\\% \$")

legend(loc = "upper center", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 4)
savefig(string("plots_final//delta_up.pdf"), format = :pdf)

##
##

using PyPlot

fig = figure(figsize=(8, 2.2))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.09,0.22,0.905,0.74])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)
xticks([10,20,30,40,50], labels = ["10", "20", "40", "50"])

#yticks([0.001,0.01,0.1,0.0,1.0,1000.0], labels = ["\&10^{-3}\&","\&10^{-2}\&","\&10^{-1}\&","\&0\&","\&10^{1}\&","\&10^{3}\&"])
#yticks([0.01, 0.0], labels = ["10", "0"])

#ax.set_yscale("log")
#ax.set_axisbelow(true)
xlabel("\$g\$")
ylim(bottom = 0.1, top = 1250)
#xlim(left = 0.0, right = length(u_buses))
ylabel("\$\\underline{\\delta}[\\\$]\$")

plot(ges, -cc2_41, color = "crimson", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~4.85\\% \$")
plot(ges, -cc2_42, color = "yellow", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~9.8\\% \$")
plot(ges, -cc2_43, color = "orange", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~19.76\\% \$")
plot(ges, -cc2_44, color = "silver", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~38.98\\% \$")

legend(loc = "upper center", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 4)
savefig(string("plots_final//delta_down.pdf"), format = :pdf)

##
##

# fig = figure(figsize=(8, 2.8))
# rc("font", family = "serif", style = "italic", size = 14)
# rc("text", usetex = true)
# rc("lines", linewidth = 1)
#
# ax = fig.add_axes([0.1,0.17,0.89,0.815])
# grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
# ax.tick_params(direction = "in", top = true, right = true, width = 1.4)
#
# #yticks([0.001,0.01,0.1,0.0,1.0,1000.0], labels = ["\&10^{-3}\&","\&10^{-2}\&","\&10^{-1}\&","\&0\&","\&10^{1}\&","\&10^{3}\&"])
# #yticks([0.01, 0.0], labels = ["10", "0"])
#
# #ax.set_yscale("log")
# #ax.set_axisbelow(true)
# xlabel("\$g\$")
# #ylim(bottom = 0.1, top = 60000000)
# #xlim(left = 0.0, right = length(u_buses))
# ylabel("\$[\\\$]\$")
# lmps44
# plot(ges, -cc1_11, color = "crimson", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~4.85\\% \$")
# plot(ges, -cc1_12, color = "yellow", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~9.8\\% \$")
# plot(ges, -cc1_13, color = "orange", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~19.76\\% \$")
# plot(ges, -cc1_14, color = "silver", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~38.98\\% \$")
# # plot(string(u_buses), [γpc4 for i in u_buses], color = "crimson", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~4.85\\% \$")
# #
# # plot(string(u_buses), χ1up_cheb3, color = "navy", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~9.7\\% \$")
# # plot(string(u_buses), [γpc3 for i in u_buses], color = "navy", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~9.7\\% \$")
# #
# # plot(string(u_buses), χ1up_cheb2, color = "silver", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~19.38\\% \$")
# # plot(string(u_buses), [γpc2 for i in u_buses], color = "silver", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~19.38\\% \$")
# #
# # plot(string(u_buses), χ1up_cheb1, color = "limegreen", lw = 1.2, ls = "dotted", marker = "D", ms = 3.0, mfc = "white", label = "\$\\chi,~38.8\\% \$")
# # plot(string(u_buses), [γpc1 for i in u_buses], color = "limegreen", lw = 1.2, ls = "dashed", mfc = "white", label = "\$\\gamma,~38.8\\% \$")
#
# legend(loc = "upper center", fancybox = false, edgecolor = "black", framealpha = 0.9, ncol = 4)
# savefig(string("plots_final//lmps3.pdf"), format = :pdf)
