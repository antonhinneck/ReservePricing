using PyPlot

u_buses = Vector{String}()
s_lmp_sys = Vector{Float64}()
s_lmp_n2n = Vector{Float64}()
σ = Vector{Float64}()

for farm in farms

    push!(u_buses, string(farm.bus))
    push!(σ, farm.σ)
    push!(s_lmp_sys, λ_sys[farm.bus])
    push!(s_lmp_n2n, λ_n2n[farm.bus])

end

cd(@__DIR__)
fig = figure(figsize=(8, 3.2))
rc("font",family="serif",style="italic")
rc("mathtext",fontset="dejavuserif")
rc("lines",linewidth=1)

ax = fig.add_axes([0.09,0.125,0.9,0.86])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
#ax.tick_params(direction="in",top=true,right=true,width=1.4)

#ax.set_axisbelow(true)
xlabel("\$u\$")
ylabel("\$\\lambda_{u}\$")
##xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

plot(u_buses, s_lmp_sys, color = "blue",  label = "node-to-node", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)
plot(u_buses, s_lmp_n2n, color = "orange",  label = "node-to-node", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)

savefig(string("lmp_sys_n2n.pdf"), format = :pdf)

using PyPlot
cd(@__DIR__)
fig = figure(figsize=(8, 2.6))
rc("font",family="serif",style="italic")
rc("mathtext",fontset="dejavuserif")
rc("lines",linewidth=1)

ax = fig.add_axes([0.08,0.16,0.91,0.82])
grid(linewidth = 0.2, linestyle = (0, (10, 10)), color = "lightgray")
#ax.tick_params(direction="in",top=true,right=true,width=1.4)

#ax.set_axisbelow(true)
xlabel("\$u\$")
ylabel("\$\\sigma_{u}\$")
##xlim(left=-5,right=5)
#x = [0.01 * i for i in -50000:50000]

plot(u_buses, σ, color = "navy",  label = "node-to-node", lw = 1, ls = "dashed", marker = "+", ms = 7.4, mew = 1.6)

savefig(string("variances.pdf"), format = :pdf)
