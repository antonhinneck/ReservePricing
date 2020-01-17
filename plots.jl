using PyPlot
cd(@__DIR__)

u_buses = Vector{String}()
s_lmp = Vector{Float64}()
s_lmp_n2n = Vector{Float64}()
a_lmp = Vector{Float64}()
a_lmp_n2n = Vector{Float64}()
σ = Vector{Float64}()
bus = [i for i in 1:n_buses]
gens = [string(i) for i in 1:n_generators]
ures = [i for i in 1:n_farms]

for farm in farms

    push!(u_buses, string(farm.bus))
    push!(σ, farm.σ)
    push!(s_lmp, λ[farm.bus] / 100)
    push!(s_lmp_n2n, λ_ab[farm.bus] / 100)
    push!(a_lmp, λ_n2n[farm.bus] / 100)
    push!(a_lmp_n2n, λ_n2n_ab[farm.bus] / 100)

end

fig = figure(figsize=(8, 3))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.09,0.14,0.9,0.85])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction="in",top=true,right=true,width=1.4)

#ax.set_axisbelow(true)
xlabel("\$u\$")
ylabel("\$\\lambda_{u}\$")
##xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

plot(u_buses, s_lmp, color = "lightblue", mec = "blue", mfc = "blue", label = "sym, system-wide", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)
plot(u_buses, s_lmp_n2n, color = "lightgreen", mec = "green", mfc = "green", label = "sym, node-to-node", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)
plot(u_buses, a_lmp, color = "yellow", mec = "orange", mfc = "orange", label = "asym, system-wide", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)
plot(u_buses, a_lmp_n2n, color = "coral", mec = "red", mfc = "red", label = "asym, node-to-node", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)

legend(loc = "lower right",fancybox=false, edgecolor="black")
savefig(string("plots//lmps_u.pdf"), format = :pdf)

fig = figure(figsize=(8, 3))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

## SYSTEM-WIDE AND N2N
##--------------------

ax = fig.add_axes([0.09,0.14,0.9,0.85])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction="in",top=true,right=true,width=1.4)

#ax.set_axisbelow(true)
xlabel("\$i\$")
ylabel("\$\\lambda_{i}\$")
##xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

plot(bus, λ / 100, color = "lightblue", mec = "blue", mfc = "blue", label = "sym, system-wide", lw = 0.8, ls = "dashed", marker = "D", ms = 1.5, mew = 1)
plot(bus, λ_n2n / 100, color = "lightgreen", mec = "green", mfc = "green", label = "sym, node-to-node", lw = 0.8, ls = "dashed", marker = "D", ms = 1.5, mew = 1)
plot(bus, λ_ab / 100, color = "yellow", mec = "orange", mfc = "orange", label = "asym, system-wide", lw = 0.8, ls = "dashed", marker = "D", ms = 1.5, mew = 1)
plot(bus, λ_n2n_ab / 100, color = "coral", mec = "red", mfc = "red", label = "asym, node-to-node", lw = 0.8, ls = "dashed", marker = "D", ms = 1.5, mew = 1)

legend(loc = "lower center",fancybox=false, edgecolor="black")
savefig(string("plots//lmps.pdf"), format = :pdf)

  # VARIANCES

fig = figure(figsize=(8, 2.6))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.09,0.18,0.9,0.8])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction="in",top=true,right=true,width=1.4)

#ax.set_axisbelow(true)
xlabel("\$u\$")
ylabel("\$\\sigma_{u}\$")
ylim(bottom=0.05,top=0.45)
#x = [0.01 * i for i in -50000:50000]

plot(u_buses, σ, color = "silver", mec = "red", mfc = "red",  label = "\$\\sigma_{u}\$", lw = 1, ls = "dashed", marker = "D", ms = 4, mew = 1)

legend(loc = "upper right",fancybox=false, edgecolor="black")
savefig(string("plots//variances.pdf"), format = :pdf)

## SYSTEM-WIDE
##------------

fig = figure(figsize=(5, 3))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.14,0.15,0.82,0.84])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction="in",top=true,right=true,width=1.4)

#ax.set_axisbelow(true)
xlabel("Case")
ylabel("\$\\gamma\$")
##xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

plot(["sym", "asym -", "asym +"], [γ, γm, γp], color = "lightblue", mec = "blue", mfc = "blue", label = "symmetric", lw = 0.8, ls = "dashed", marker = "D", ms = 2.5, mew = 1)

savefig(string("plots//gamma.pdf"), format = :pdf)

## N2N
##----

fig = figure(figsize=(8, 3))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.09,0.15,0.9,0.84])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction = "in", top = true, right = true, width = 1.4)

#ax.set_axisbelow(true)
xlabel("\$u\$")
ylabel("\$\\chi^{(\\pm)}_{u}\$")
##xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

plot(u_buses, χm, color = "teal", mec = "navy", mfc = "navy", label = "Model 4, \$\\chi^{-}_{u}\$", lw = 1, ls = "dashed", marker = "D", ms = 3, mew = 1)
plot(u_buses, χp, color = "lightgreen", mec = "blue", mfc = "blue", label = "Model 4, \$\\chi^{+}_{u}\$", lw = 1, ls = "dashed", marker = "D", ms = 3, mew = 1)
plot(u_buses, χ, color = "gray", mec = "black", mfc = "black", label = "Model 2, \$\\chi_{u}\$", lw = 1, ls = "dashed", marker = "D", ms = 3, mew = 1)

legend(loc = "lower left", fancybox = false, edgecolor = "black")
savefig(string("plots//chi.pdf"), format = :pdf)

## Participation Factors
##----------------------

fig = figure(figsize=(8, 3))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.09,0.18,0.9,0.8])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction="in",top=true,right=true,width=1.4)

#ax.set_axisbelow(true)
xlabel("\$i\$")
ylabel("\$\\alpha_{i}\$")
#ylim(bottom=0.05,top=0.45)
#x = [0.01 * i for i in -50000:50000]

plot(gens, a_s, color = "lightgray", mec = "navy", mfc = "navy",  label = "\$\\sigma_{u}\$", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)
plot(gens, ap, color = "yellow", mec = "red", mfc = "red",  label = "\$\\sigma_{u}\$", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)
plot(gens, am, color = "yellow", mec = "red", mfc = "red",  label = "\$\\sigma_{u}\$", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)
plot(gens, a_n2n, color = "lightblue", mec = "blue", mfc = "blue",  label = "\$\\sigma_{u}\$", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)
plot(gens, ap_n2n, color = "lightgreen", mec = "green", mfc = "green",  label = "\$\\sigma_{u}\$", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)
plot(gens, am_n2n, color = "lightgreen", mec = "green", mfc = "green",  label = "\$\\sigma_{u}\$", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)

legend(loc = "upper right",fancybox=false, edgecolor="black")
savefig(string("plots//alphas.pdf"), format = :pdf)

## Costs
##------

fig = figure(figsize=(8, 3))
rc("font", family = "serif", style = "italic", size = 14)
rc("text", usetex = true)
rc("lines", linewidth = 1)

ax = fig.add_axes([0.09,0.18,0.9,0.8])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
ax.tick_params(direction="in",top=true,right=true,width=1.4)

#ax.set_axisbelow(true)
xlabel("\$i\$")
ylabel("\$c_{i}\$")
#ylim(bottom=0.05,top=0.45)
#x = [0.01 * i for i in -50000:50000]

plot(gens, c_vec, color = "lightgray", mec = "navy", mfc = "navy",  label = "\$\\sigma_{u}\$", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)

legend(loc = "upper right",fancybox=false, edgecolor="black")
savefig(string("plots//costs.pdf"), format = :pdf)
