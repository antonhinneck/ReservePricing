mutable struct Generator
   index::Any
   bus_idx::Int
   g_max::Float64
   cost::Float64
   quad_cost::Float64
   function Generator(index, bus_idx, g_max)
      g = new()
      g.index  = index
      g.bus_idx = bus_idx
      g.cost = 1.00
      g.quad_cost = 0.00
      g.g_max = g_max
      return g
   end
end

mutable struct Bus
   index::Any
   is_slack::Bool
   d_P::Float64
   d_Q::Float64
   gen_list::Vector{Int}
   farm_list::Vector{Int}
   start_line::Vector{Int}
   end_line::Vector{Int}
   cosphi::Float64
   tanphi::Float64
   function Bus(index, d_P, d_Q, is_slack)
      b = new()
      b.index = index
      b.is_slack = is_slack
      b.d_P = d_P
      b.d_Q = d_Q
      b.gen_list = []
      b.farm_list = []
      cosphi = d_P/(sqrt(d_P^2 + d_Q^2))
      tanphi = d_Q/d_P
      if isnan(cosphi)
        b.cosphi = 0
        b.tanphi = 0
      else
        b.cosphi = cosphi
        b.tanphi = tan(acos(cosphi))
      end
      b.start_line = []
      b.end_line = []
      return b
   end
end
function set_bus_active_load(b::Bus, dP; auto_set_dQ=false)
    b.d_P = dP
    if auto_set_dQ
        b.d_Q = dP * b.tanphi
    end
end

mutable struct Line
   index::Any
   to_node::Int # the "to" node
   from_node::Int # the "from" node
   r::Float64 # the resistance value
   x::Float64 # the reactance value
   b::Float64 # the susceptance value
   s_max::Float64 # the capacity of the line
   function Line(index, to_node, from_node, r, x, s_max; b=-1)
      l = new()
      l.index = index
      l.to_node = to_node
      l.from_node = from_node
      l.r = r
      l.x = x
      if b == -1
        l.b = (x/(r^2 + x^2))
      else
        l.b = b
      end
      l.s_max = s_max
      return l
   end
end


mutable struct Farm
    μ::Float64
    σ::Float64
    bus::Int
end


function load_network(datadir)

    nodes_raw = CSV.read("$datadir/nodes.csv")
    length(unique(nodes_raw[:index]))==length(nodes_raw[:index]) ? nothing : @warn("Ambiguous Node Indices")

    lines_raw = CSV.read("$datadir/lines.csv")
    length(unique(lines_raw[:index]))==length(lines_raw[:index]) ? nothing : @warn("Ambiguous Line Indices")

    generators_raw = CSV.read("$datadir/plants.csv")
    length(unique(generators_raw[:index]))==length(generators_raw[:index]) ? nothing : @warn("Ambiguous Generator Indices")

    println(">>> Loading Buses")
    buses = []
    for n in 1:size(nodes_raw,1)
        index = nodes_raw[n, :index]
        d_P = nodes_raw[n, :Pd]
        d_Q = nodes_raw[n, :Qd]
        is_slack = nodes_raw[n, :slack] == "True" ? true : false
        newb = Bus(index, d_P, d_Q, is_slack)
        push!(buses, newb)
    end

    println(">>> Loading Lines")
    lines = []
    for l in 1:size(lines_raw, 1)
        index = lines_raw[l, :index]
        from_node = lines_raw[l, :node_i]
        to_node = lines_raw[l, :node_j]
        r = lines_raw[l, :r]
        x = lines_raw[l, :x]
        b = lines_raw[l, :b]
        s_max = lines_raw[l, :maxflow]
        newl = Line(index, to_node, from_node, r, x, s_max; b=b)
        push!(buses[newl.from_node].start_line, index)
        push!(buses[newl.to_node].end_line, index)
        push!(lines, newl)
    end

    println(">>> Loading Generators")
    generators = []
    for g in 1:size(generators_raw, 1)
        index = generators_raw[g, :index]
        bus_idx = generators_raw[g, :node]
        g_max = generators_raw[g, :g_max]
        cost = generators_raw[g, :mc]
        newg = Generator(index, bus_idx, g_max)
        newg.cost = cost
        push!(buses[newg.bus_idx].gen_list, newg.index)
        push!(generators, newg)
    end

    println(">>> Loaded $(size(buses, 1)) buses, $(size(lines, 1)) lines and $(size(generators, 1)) generators.")

    return buses, lines, generators

end
