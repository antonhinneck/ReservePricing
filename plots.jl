using PyPlot
cd(@__DIR__)

u_buses = Vector{String}()
s_lmp_sys = Vector{Float64}()
s_lmp_n2n = Vector{Float64}()
σ = Vector{Float64}()

for farm in farms

    push!(u_buses, string(farm.bus))
    push!(σ, farm.σ)
    push!(s_lmp_sys, λ_sys[farm.bus] / 100)
    push!(s_lmp_n2n, λ_n2n[farm.bus] / 100)

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

plot(u_buses, s_lmp_sys, color = "lightblue", mec = "blue", mfc = "blue", label = "system-wide", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)
plot(u_buses, s_lmp_n2n, color = "lightgreen", mec = "green", mfc = "green", label = "node-to-node", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)

legend(loc = "lower right",fancybox=false, edgecolor="black")
savefig(string("plots//lmp_sys_n2n.pdf"), format = :pdf)

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

plot(u_buses, σ, color = "lightgray", mec = "navy", mfc = "navy",  label = "\$\\sigma_{u}\$", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)

legend(loc = "upper right",fancybox=false, edgecolor="black")
savefig(string("plots//variances.pdf"), format = :pdf)
