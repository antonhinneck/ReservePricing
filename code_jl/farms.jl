#wind_buses = [3, 8, 11, 20, 24, 26, 31, 38, 43, 49, 53]
#wind_cpcty = [70.0, 147.0, 102.0, 105.0, 113.0, 84.0, 59.0, 250.0, 118.0, 76.0, 72.0]

function create_wind_farms(;buses = [3, 8, 11, 20, 24, 26, 31, 38, 43, 49, 53], capacity = [70.0, 147.0, 102.0, 105.0, 113.0, 84.0, 59.0, 250.0, 118.0, 76.0, 72.0], scaling_sigma = 1.0, scaling_cap = 1.0)

    @assert length(buses) == length(capacity)
    farms = Vector{Farm}()
    nf = length(buses)
    capacity = capacity * scaling_cap

    for i in 1:nf
        push!(farms,  Farm(capacity[i] / 100, scaling_sigma * capacity[i] / 10 / 100, buses[i]))
    end

    σ_vec = [i.σ^2 for i in farms]

    Σ = diagm(0 => (σ_vec))
    s_sq = sum(Σ)
    Σ_rt = sqrt(Σ)
    s = sum(Σ_rt)

    return farms, nf, σ_vec, Σ, s_sq, Σ_rt, s
end
