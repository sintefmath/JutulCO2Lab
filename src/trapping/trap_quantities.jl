# Quantities derived from a trap analysis: trap volumes, downstream traps,
# capacity trees, the best injection point and the flattened caprock
# (co2lab's volumesOfTraps, downstreamTraps, maximizeTrapping,
# findOptimalInjectionPoint, flattenTraps).

# ---------------------------------- Volumes -----------------------------------

"""
    trap_volumes(s, res; traps = 1:ntrap, poro = ones) -> Vector{Float64}

Pore volume of each trap (co2lab's `volumesOfTraps`): the space between the
caprock and the trap's spill depth, cut at the formation thickness.

# Arguments

- `s`: the [`trap_surface`](@ref).
- `res`: its [`trap_analysis`](@ref).
- `traps`: the traps to measure; default all.
- `poro`: porosity per column; default 1 (bulk volume).

# Returns

A new vector, one volume (m³) per trap in `traps`:
Σ area · min(spill depth − z, H) · poro over the trap's columns.

# Notes

- co2lab also adds `eps` m³ per column; this does not.
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
`downstreamTraps`), from the [`trap_analysis`](@ref) result `res`.

# Arguments

- `res`: the [`trap_analysis`](@ref) result.
- `t`: the trap number.

# Returns

A new vector of trap numbers, starting with `t`, then ordered by the number
of spill steps from `t`, and within a step by increasing trap number.
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

The spill trees with the largest trapping capacity (co2lab's
`maximizeTrapping`): where injected CO2 could fill the most trap volume.

A tree is a root trap and everything downstream of it; its value is the sum
of their [`trap_volumes`](@ref) (unit porosity). Trees are picked greedily
by value.

# Arguments

- `s`: the [`trap_surface`](@ref).
- `res`: its [`trap_analysis`](@ref).
- `n`: the maximum number of trees; default all.
- `remove_overlap`: count a trap only in the first tree that contains it.
- `calculate_all`: use every trap as a root; default only the traps nothing
  spills into.

# Returns

A new vector of NamedTuples, largest value first:
- `root`: the root trap;
- `traps`: the trap numbers of the tree, sorted;
- `value`: its volume (m³), less traps counted earlier if `remove_overlap`.

# Notes

- Equal values are taken in root order, as co2lab's `max`.
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

The column where injected CO2 can fill the most trap volume (co2lab's
`findOptimalInjectionPoint`). Candidates are the interfaces between two
spill regions, or between a spill region and the outside; CO2 injected on
one reaches the traps downstream of both sides.

# Arguments

- `s`: the [`trap_surface`](@ref).
- `res`: its [`trap_analysis`](@ref).

# Returns

A new NamedTuple:
- `column`: the shallowest column on the best interface;
- `volume`: the trap volume (m³, unit porosity) reachable from it;
- `edges`: the surface edges of the interface (indices into `s.edges`);
- `traps`: the traps of its two regions (0: outside).

# Notes

- Sealed edges are not interfaces to the outside.
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
`flattenTraps`): trap columns and their nodes are moved to the spill depth.

# Arguments

- `s`: the [`trap_surface`](@ref).
- `res`: its [`trap_analysis`](@ref).

# Returns

A new surface NamedTuple, `s` with new `zc` and `z`; the other fields are
shared with `s`.

# Notes

- A node shared by trap columns of different depth takes the last column's
  depth, as in co2lab.
"""
function flatten_traps(s, res)
    zc, z = copy(s.zc), copy(s.z)
    for c in findall(>(0), res.traps)
        zc[c] = res.trap_z[res.traps[c]]
        z[s.polygons[c]] .= zc[c]
    end
    return (; s..., zc = zc, z = z)
end
