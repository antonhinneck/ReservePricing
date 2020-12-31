#wind_buses = [3, 8, 11, 20, 24, 26, 31, 38, 43, 49, 53]
#wind_cpcty = [70.0, 147.0, 102.0, 105.0, 113.0, 84.0, 59.0, 250.0, 118.0, 76.0, 72.0]

mutable struct uRES

    forecast::Float64
    μ::Float64
    σ::Float64
    bus::Int64

end

function create_wind_farms(; scaling_sigma = 1.0, scaling_cap = 1.0, fc = 1)

    capacity = 2 .* [70.0, 147.0, 102.0, 105.0, 113.0, 84.0, 59.0, 250.0, 118.0, 76.0, 72.0] ./ 100
    bus_ids = [3, 8, 11, 20, 24, 26, 31, 38, 43, 49, 53]

    @assert length(bus_ids) == length(capacity)
    uRESs = Vector{uRES}()
    nf = length(bus_ids)
    capacity = capacity

    ud = JLD.load("uncertainty_datafile.jld")["data"]

    for i in 1:nf
        push!(uRESs,  uRES(capacity[i] * ud["fcs$(fc)"][i], capacity[i] * ud["means"][i], capacity[i] * ud["stdDevs"][i], bus_ids[i]))
        global buses[bus_ids[i]].farmids = Vector{Int64}()
        push!(buses[bus_ids[i]].farmids, i)
    end

    σ = [i.σ for i in uRESs]
    μ = [i.μ for i in uRESs]

    Σ = diagm(0 => (σ.^2))
    s_sq = sum(Σ)
    Σ_rt = sqrt(Σ)
    s = sqrt(s_sq)
    println(string("s: ",s,"---s: ",sum(Σ_rt)))

    return uRESs, nf, σ, Σ, s_sq, Σ_rt, s, μ
end

# uncertainty_data = JLD.load("uncertainty_datafile.jld")
# uncertainty_data["data"]["means"]
# uncertainty_data["data"]["stdDevs"]
# uncertainty_data["data"]["fcs5"]
