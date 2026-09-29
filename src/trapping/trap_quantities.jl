# Quantities derived from a trap analysis (roadmap Stage 8 item 5): trap
# volumes, downstream traps, capacity trees, the best injection point and the
# flattened caprock — co2lab's volumesOfTraps, downstreamTraps,
# maximizeTrapping, findOptimalInjectionPoint, flattenTraps.

# ---------------------------------- Volumes -----------------------------------

"""
    trap_volumes(s, res; traps = 1:ntrap, poro = ones) -> Vector{Float64}

Pore volume of each trap in `traps`: Σ area · min(spill depth − z, H) · poro
over its columns (co2lab's `volumesOfTraps`, which also adds `eps` m³ per
column; not here). `s` is the [`trap_surface`](@ref), `res` the
[`trap_analysis`](@ref) result, `poro` per column. Returns a new vector.
"""
function trap_volumes(s, res; traps = 1:maximum(res.traps; init = 0),
                      poro = ones(length(s.zc)))
    return map(traps) do t
        cells = findall(==(t), res.traps)
        # CO2 column from the spill depth up, cut at the formation thickness
        h = min.(res.trap_z[t] .- s.zc[cells], s.H[cells])
        sum(s.area[cells] .* h .* poro[cells])
    end
end

# ------------------------------ Downstream traps ------------------------------

"""
    downstream_traps(res, t) -> Vector{Int}

Trap `t` and every trap it spills into, directly or further on (co2lab's
`downstreamTraps`): in order of the number of spill steps, each step in
increasing trap number. Returns a new vector.
"""
function downstream_traps(res, t::Integer)
    A = res.trap_adj
    out, frontier = [t], [t]
    # walks of length 1, 2, … (the spill graph has no cycles)
    for _ in 1:size(A, 1)
        frontier = sort!(unique([j for i in frontier for j in findall(A[i, :])]))
        isempty(frontier) && break
        append!(out, frontier)
    end
    return unique(out)
end

# ------------------------------- Capacity trees -------------------------------

"""
    maximize_trapping(s, res; n = Inf, remove_overlap = true,
                      calculate_all = false) -> Vector{NamedTuple}

co2lab's `maximizeTrapping`: the spill trees with the largest trapping
capacity. A tree is a root trap (by default the traps nothing spills into;
all traps with `calculate_all`) and everything downstream of it; its value
is the sum of their [`trap_volumes`](@ref) (unit porosity). Trees are picked
greedily by value, at most `n`; with `remove_overlap` a trap counts only in
the first tree that contains it. Returns a new vector of `(root, traps,
value)`, largest first.
"""
function maximize_trapping(s, res; n = Inf, remove_overlap::Bool = true,
                           calculate_all::Bool = false)
    v = trap_volumes(s, res)
    nt = length(v)
    # start points: every trap, or only those nothing spills into
    roots = calculate_all ? collect(1:nt) :
            [t for t in 1:nt if !any(res.trap_adj[:, t])]
    down = [sort(downstream_traps(res, r)) for r in roots]
    trees = NamedTuple{(:root, :traps, :value), Tuple{Int, Vector{Int}, Float64}}[]
    done = fill(false, length(roots))
    for _ in 1:min(n, length(roots))
        vals = [done[k] ? 0.0 : sum(v[down[k]]) for k in eachindex(roots)]
        # first of equal maxima, as co2lab's max
        k = argmax(vals)
        push!(trees, (root = roots[k], traps = down[k], value = vals[k]))
        # a trap filled by this tree counts no more for later ones
        remove_overlap && (v[down[k]] .= 0.0)
        done[k] = true
    end
    return trees
end

# ------------------------------ Injection point -------------------------------

"""
    optimal_injection_point(s, res) -> NamedTuple

co2lab's `findOptimalInjectionPoint`: among the interfaces between two spill
regions (or a spill region and the outside), the one whose two sides can
fill the most trap volume downstream, and on it the shallowest column.
Returns `(column, volume, edges, traps)`: the column, that volume, the
surface edges of the interface and the two regions' traps (0: outside).
"""
function optimal_injection_point(s, res)
    v = trap_volumes(s, res)
    region(c) = c == 0 ? 0 : res.trap_regions[c]
    # a sealed edge is no interface to the outside
    pairs = [_is_sealed(s, k) ? (-1, -1) :
             Tuple(sort([region(s.edge_columns[1, k]), region(s.edge_columns[2, k])]))
             for k in axes(s.edge_columns, 2)]
    # candidates: every pair of different regions that share an edge
    candidates = sort!(unique(filter(p -> p[1] != p[2], pairs)))
    reach(p) = unique(vcat([downstream_traps(res, t) for t in p if t > 0]...))
    volumes = [sum(v[reach(p)]) for p in candidates]
    best = candidates[argmax(volumes)]
    edges = findall(==(best), pairs)
    cols = filter(>(0), sort!(unique(vec(s.edge_columns[:, edges]))))
    # the shallowest column on the best interface
    return (column = cols[argmin(s.zc[cols])], volume = maximum(volumes),
            edges = edges, traps = best)
end

# ---------------------------------- Flatten -----------------------------------

"""
    flatten_traps(s, res) -> NamedTuple

The surface with every trap filled to its spill depth (co2lab's
`flattenTraps`): trap columns and their nodes at the spill depth (a node
shared by columns of different depth takes the last column's, as co2lab).
Returns a new surface; `s` is only read.
"""
function flatten_traps(s, res)
    zc, z = copy(s.zc), copy(s.z)
    for c in findall(>(0), res.traps)
        zc[c] = res.trap_z[res.traps[c]]
        z[s.polygons[c]] .= zc[c]
    end
    return (; s..., zc = zc, z = z)
end
