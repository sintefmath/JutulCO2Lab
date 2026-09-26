# Structural traps from column depths and column adjacency: co2lab-spillpoint's
# cell-based method, trapAnalysis(Gt, true) — findTrappingStructure,
# findTrapConnections, maxTPFAGravityMatrix (roadmap Stage 9 item 2).
#
# Depths are positive down: CO2 rises to smaller z. A "top" column has no
# shallower neighbour; a trap is the set of columns reached downhill from a top
# that lie above the shallowest column on the rim of the reached set (the
# spill point). Traps are filled to their spill depth level by level, so that
# traps spilling into each other merge.

# --------------------------------- Adjacency ----------------------------------

# Neighbour lists of the columns (sorted), and whether each column touches the
# domain boundary, from the surface's edges
function _column_neighbours(s)
    n = length(s.zc)
    nbrs = [Int[] for _ in 1:n]
    on_boundary = fill(false, n)
    for k in axes(s.edge_columns, 2)
        a, b = s.edge_columns[1, k], s.edge_columns[2, k]
        if a > 0 && b > 0
            push!(nbrs[a], b); push!(nbrs[b], a)
        else
            on_boundary[max(a, b)] = true
        end
    end
    foreach(sort!, nbrs)
    return (nbrs, on_boundary)
end

# Columns with at least one neighbour and none shallower (co2lab's "top")
_top_columns(z, nbrs) =
    [c for c in eachindex(z) if !isempty(nbrs[c]) && all(z[b] >= z[c] for b in nbrs[c])]

# Columns reached from `start` going down or level (z not decreasing)
function _downhill_reach(nbrs, z, start)
    reached = fill(false, length(z))
    reached[start] = true
    stack = [start]
    while !isempty(stack)
        c = pop!(stack)
        for b in nbrs[c]
            if !reached[b] && z[b] >= z[c]
                reached[b] = true
                push!(stack, b)
            end
        end
    end
    return reached
end

# ---------------------------- Trapping structure ------------------------------

# One top's trap in the current (partly filled) depths `zf` (co2lab's
# bfs_argumented): the reached columns at or above the shallowest rim column,
# and that spill depth; no trap if none lies strictly above it
function _trap_from_top(nbrs, on_boundary, zf, zlevel, top)
    reached = _downhill_reach(nbrs, zlevel, top)
    # rim: reached columns next to an unreached one or the domain boundary
    rim = [c for c in findall(reached)
           if on_boundary[c] || any(b -> !reached[b], nbrs[c])]
    z_spill = minimum(zf[rim])
    cells = [c for c in findall(reached) if zf[c] <= z_spill]
    any(c -> zf[c] < z_spill, cells) || return (Int[], z_spill)
    return (cells, z_spill)
end

# co2lab's findTrappingStructure: the tops, each column's spill depth (0 if
# not in a trap) and the filled depths
function _trapping_structure(s, nbrs, on_boundary)
    z = s.zc
    tops = _top_columns(z, nbrs)
    zf = copy(z)
    candidate = fill(false, length(z)); candidate[tops] .= true
    z_spill_loc = zeros(length(z))
    while any(candidate)
        # the downhill graph is fixed per level; the depths fill as traps are found
        zlevel = copy(zf)
        level = Tuple{Vector{Int}, Float64}[]
        for t in tops
            candidate[t] || continue
            cells, z_spill = _trap_from_top(nbrs, on_boundary, zf, zlevel, t)
            if isempty(cells)
                candidate[t] = false
            else
                candidate[cells] .= false; candidate[t] = true
                zf[cells] .= z_spill
                push!(level, (cells, z_spill))
            end
        end
        for (cells, z_spill) in level
            z_spill_loc[cells] .= max.(z_spill_loc[cells], z_spill)
        end
    end
    # columns exactly at their spill depth are not in the trap
    z_spill_loc[z_spill_loc .== z] .= 0.0
    return (tops = filter(t -> z_spill_loc[t] > 0, tops), z_spill_loc = z_spill_loc,
            zf = zf)
end

# ------------------------------ Traps and rivers ------------------------------

# Trap index per column (0: none): columns with the same spill depth, split
# into connected pieces (co2lab's processPartition; pieces numbered by their
# first column here, where co2lab's dmperm order is arbitrary)
function _number_traps(nbrs, z_spill_loc)
    traps = zeros(Int, length(z_spill_loc))
    n = 0
    for c in eachindex(z_spill_loc)
        (z_spill_loc[c] > 0 && traps[c] == 0) || continue
        n += 1
        traps[c] = n
        stack = [c]
        while !isempty(stack)
            a = pop!(stack)
            for b in nbrs[a]
                if traps[b] == 0 && z_spill_loc[b] == z_spill_loc[c]
                    traps[b] = n
                    push!(stack, b)
                end
            end
        end
    end
    return traps
end

# co2lab's findTrapConnections' rivers: from each trap (deepest spill first),
# out over its shallowest rim neighbour and up the steepest way until another
# trap, the boundary or a local top. Returns (from, to, cells) per river, `to`
# 0 when it leaves the domain or ends at a top.
function _trace_cell_rivers(nbrs, on_boundary, z, traps)
    ntrap = maximum(traps; init = 0)
    trap_z = [z[findfirst(==(t), traps)] for t in 1:ntrap]
    visited = fill(false, ntrap)
    rivers = Tuple{Int, Int, Vector{Int}}[]
    while !all(visited)
        # deepest unvisited trap first (co2lab: max of traps_z)
        t = argmax(i -> visited[i] ? -Inf : trap_z[i], 1:ntrap)
        while t > 0 && !visited[t]
            visited[t] = true
            r = _river_from_trap(nbrs, on_boundary, z, traps, t)
            isnothing(r) && break
            push!(rivers, (t, r.to, r.cells))
            t = r.to
        end
    end
    return rivers
end

# One river out of trap t, or nothing if no rim neighbour is shallower than
# or level with the spill depth
function _river_from_trap(nbrs, on_boundary, z, traps, t)
    cells = findall(==(t), traps)
    rim = sort!(unique([b for c in cells for b in nbrs[c] if traps[b] != t]))
    up = filter(b -> z[b] <= z[cells[1]], rim)
    isempty(up) && return nothing
    cell = up[argmin(z[up])]
    line = [cell]
    while !on_boundary[cell]
        cc = sort!([cell; nbrs[cell]])
        next = cc[argmin(z[cc])]
        z[next] >= z[cell] && break
        cell = next
        traps[cell] > 0 && break
        push!(line, cell)
    end
    return (to = traps[cell], cells = line)
end

# ------------------------------- Spill regions --------------------------------

# Each column's steepest-ascent neighbour (0 if none is shallower), and the
# level neighbours (co2lab's maxTPFAGravityMatrix: slope over the distance
# between centroids; ties to the lowest column index)
function _steepest_ascent(s, nbrs)
    z, up = s.zc, zeros(Int, length(s.zc))
    for c in eachindex(z)
        best = 0.0
        for b in nbrs[c]
            # co2lab's distance, sqrt(dx² + dy²) (hypot may differ in the last bit)
            slope = (z[c] - z[b]) / sqrt((s.cx[c] - s.cx[b])^2 + (s.cy[c] - s.cy[b])^2)
            slope > best && ((best, up[c]) = (slope, b))
        end
    end
    return up
end

# Columns draining to `top`: those whose steepest-ascent path, moving freely
# between level neighbours, reaches it (co2lab: bfs on the transposed matrix)
function _drainage(nbrs, z, up, top)
    down = [Int[] for _ in eachindex(z)]
    for c in eachindex(up)
        up[c] > 0 && push!(down[up[c]], c)
    end
    reached = fill(false, length(z)); reached[top] = true
    stack = [top]
    while !isempty(stack)
        c = pop!(stack)
        for b in [down[c]; filter(b -> z[b] == z[c], nbrs[c])]
            reached[b] || (reached[b] = true; push!(stack, b))
        end
    end
    return reached
end

# --------------------------------- Analysis -----------------------------------

"""
    _cell_trap_analysis(s) -> NamedTuple

co2lab's `trapAnalysis(Gt, true)` on a [`trap_surface`](@ref) `s`: traps
(trap index per column, 0 outside traps), their spill depths `trap_z`, the
spill region of each column (`trap_regions`: the trap it drains into, 0 if
none), which trap spills into which (`trap_adj`, a sparse `Bool` matrix),
the rivers out of each trap (`rivers[t]`: vectors of columns) and the top
columns of the traps (`top`). `s` is only read. Returns a new NamedTuple.
"""
function _cell_trap_analysis(s)
    nbrs, on_boundary = _column_neighbours(s)
    st = _trapping_structure(s, nbrs, on_boundary)
    traps = _number_traps(nbrs, st.z_spill_loc)
    ntrap = maximum(traps; init = 0)
    trap_z = [st.z_spill_loc[findfirst(==(t), traps)] for t in 1:ntrap]

    # rivers are traced on the filled surface, traps at their spill depth
    zr = copy(st.zf)
    zr[st.z_spill_loc .> 0] .= st.z_spill_loc[st.z_spill_loc .> 0]
    river_list = _trace_cell_rivers(nbrs, on_boundary, zr, traps)
    rivers = [Vector{Int}[] for _ in 1:ntrap]
    foreach(r -> push!(rivers[r[1]], r[3]), river_list)
    links = [(r[1], r[2]) for r in river_list if r[2] > 0]
    trap_adj = sparse(first.(links), last.(links), true, ntrap, ntrap)

    # spill regions: what drains to each trap's top (later tops overwrite)
    up = _steepest_ascent(s, nbrs)
    trap_regions = zeros(Int, length(s.zc))
    for t in st.tops
        trap_regions[_drainage(nbrs, s.zc, up, t)] .= traps[t]
    end
    return (traps = traps, trap_z = trap_z, trap_regions = trap_regions,
            trap_adj = trap_adj, rivers = rivers, top = st.tops)
end

# ------------------------------- Front end (temporary) ------------------------

"""
    trap_analysis(s; method = :cell) -> NamedTuple

Structural traps of a [`trap_surface`](@ref) `s`, as co2lab's
`trapAnalysis(Gt, method == :cell)`. Returns [`_cell_trap_analysis`](@ref)'s
NamedTuple. `s` is only read.
"""
function trap_analysis(s; method::Symbol = :cell)
    method == :cell && return _cell_trap_analysis(s)
    throw(ArgumentError("method must be :cell (:node not yet implemented)"))
end
