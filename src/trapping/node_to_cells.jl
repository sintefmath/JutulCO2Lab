# The node-based traps projected to columns: co2lab-spillpoint's n2cTraps
# (nodeFieldToCellField, shortestPath, subregionConnectivityMatrix).
#
# - A column takes the value of its shallowest corner.
# - A trap whose columns all lie at or below its spill depth disappears; its
#   region and rivers pass on to the next remaining trap downstream.

# ------------------------------- Columns --------------------------------------

# Each column's shallowest corner (co2lab: the first of equal minima among
# west-south, east-south, west-north, east-north)
_column_rep_nodes(s) = [s.corners[argmin(s.z[s.corners[:, c]]), c]
                        for c in axes(s.corners, 2)]

# One top column per local top: among the columns whose shallowest corner is
# the top node, the shallowest (co2lab's n2c_sommets); true/false per column
function _column_tops(s, rep, top_nodes)
    id = zeros(Int, length(s.z))
    id[top_nodes] .= 1:length(top_nodes)
    best = Dict{Int, Int}()
    for c in eachindex(rep)
        k = id[rep[c]]
        k == 0 && continue
        b = get(best, k, 0)
        (b == 0 || s.zc[c] < s.zc[b]) && (best[k] = c)
    end
    tops = fill(false, length(rep))
    tops[collect(values(best))] .= true
    return tops
end

# The first remaining trap downstream of `traps` (0 if none): co2lab's
# next_remaining_downstream_trap
function _next_remaining(adj, traps, lost)
    down = sort!(unique([j for i in traps for j in findall(adj[i, :])]))
    isempty(down) && return 0
    rem = setdiff(down, lost)
    isempty(rem) || return minimum(rem)
    # every trap downstream is lost too: look one step further
    return _next_remaining(adj, intersect(down, lost), lost)
end

# Trap k (k > i) becomes k − 1, for each lost trap i in descending order
_reindex(v, lost) = [x - count(<(x), lost) for x in v]

# Trap i removed from the adjacency: its predecessors spill into its
# successors (co2lab: lmerge * C * rmerge)
function _bypass(adj, i)
    n = size(adj, 1)
    keep = setdiff(1:n, i)
    C = Float64.(adj)
    # L·C: a predecessor of i also spills wherever i spills
    L = sparse(1:n, 1:n, 1.0)[keep, :]
    L[:, i] = C[keep, i]
    return (L * C)[:, keep] .!= 0
end

# ------------------------------ Projection ------------------------------------

"""
    _node_to_cell_traps(s, nt) -> NamedTuple

co2lab's `n2cTraps`: the node-based traps `nt` ([`_node_traps`](@ref))
projected to the columns of the trap surface `s` ([`trap_surface`](@ref)).

Returns a new NamedTuple `(traps, trap_z, trap_regions, trap_adj, rivers,
top)` per column and trap, with the same fields and meaning as
[`trap_analysis`](@ref)'s result; traps left without columns are removed and
the rest renumbered.
"""
function _node_to_cell_traps(s, nt)
    rep = _column_rep_nodes(s)
    # the regions' top nodes (no downstream node)
    top_nodes = findall(n -> nt.down[n] == 0 && nt.trap_regions[n] != 0,
                        eachindex(s.z))
    regions = nt.trap_regions[rep]
    tops = _column_tops(s, rep, top_nodes)
    ntrap = length(nt.trap_z)
    traps = zeros(Int, length(rep))
    # a column is in a trap when in its region and above its spill depth
    for i in 1:ntrap
        traps[(regions .== i) .& (s.zc .< nt.trap_z[i])] .= i
    end

    # traps without columns: their regions to the next remaining one downstream
    lost = [i for i in ntrap:-1:1 if !any(==(i), traps)]
    for i in lost
        regions[regions .== i] .= _next_remaining(nt.trap_adj, i, lost)
    end
    regions, traps = _reindex(regions, lost), _reindex(traps, lost)
    trap_z = deleteat!(copy(nt.trap_z), sort(lost))
    top = sort(findall(tops .& (regions .> 0)), by = c -> regions[c])

    adj, rivers = _bypass_lost(s, nt, lost)
    look = _river_lookup(s)
    crivers = [[_truncate_river(_river_columns(s, look, r), traps) for r in rt]
               for rt in rivers]
    return (traps = traps, trap_z = trap_z, trap_regions = regions, trap_adj = adj,
            rivers = crivers, top = top)
end

# Lost traps removed from the adjacency, and the rivers into them joined to
# their rivers out (node lists, per remaining trap)
function _bypass_lost(s, nt, lost)
    adj = nt.trap_adj
    rivers = deepcopy(nt.rivers)
    for i in lost
        adj = _bypass(adj, i)
        fill_nodes = findall(==(i), nt.traps)
        for j in findall(nt.trap_adj[:, i]), k in eachindex(rivers[j])
            nt.trap_regions[rivers[j][k][end]] == i || continue
            for r in rivers[i]
                rivers[j][k] = _join_rivers(s, rivers[j][k], r, fill_nodes)
            end
        end
    end
    deleteat!(rivers, sort(lost))
    return (adj, rivers)
end

# ---------------------------------- Rivers ------------------------------------

# co2lab's join_rivers: r1 (into a lost trap) continued by r2 (out of it),
# through the trap's nodes `fill` by the shortest path if they do not meet
function _join_rivers(s, r1, r2, fill)
    r1[end] == r2[1] && return [r1[1:end-1]; r2]
    in_fill = Set(fill)
    r1 = r1[1:findlast(in(in_fill), r1)]
    r2 = r2[findfirst(in(in_fill), r2):end]
    path = _shortest_path(_fill_weights(s, fill), findfirst(==(r1[end]), fill),
                          findfirst(==(r2[1]), fill))
    return [r1[1:end-1]; fill[path[1:end-1]]; r2]
end

# co2lab's subregionConnectivityMatrix with diagonals: nodes of `fill` joined
# when they share a column, weighted by their map distance summed over the
# columns they share. Returns neighbour lists (sorted) and weights, indexed
# by position in `fill`.
function _fill_weights(s, fill)
    pos = Dict(n => k for (k, n) in enumerate(fill))
    w = Dict{Tuple{Int, Int}, Float64}()
    for p in s.polygons
        q = sort!([pos[n] for n in unique(p) if haskey(pos, n)])
        for a in eachindex(q), b in a+1:length(q)
            i, j = fill[q[a]], fill[q[b]]
            d = sqrt((s.x[i] - s.x[j])^2 + (s.y[i] - s.y[j])^2)
            w[(q[a], q[b])] = get(w, (q[a], q[b]), 0.0) + d
        end
    end
    nbrs = [Int[] for _ in fill]
    for (a, b) in keys(w)
        push!(nbrs[a], b); push!(nbrs[b], a)
    end
    foreach(sort!, nbrs)
    return (nbrs = nbrs, w = w)
end

# co2lab's shortestPath: in rounds, every node not yet reached gets the best
# distance through the reached ones (the first of equal minima); the path
# from start to target as positions
function _shortest_path(g, start, target)
    D = fill(Inf, length(g.nbrs)); D[start] = 0.0
    pred = zeros(Int, length(D))
    while isinf(D[target])
        new = Tuple{Int, Float64, Int}[]
        for u in findall(isinf, D)
            best, arg = Inf, 0
            for f in g.nbrs[u]
                isfinite(D[f]) || continue
                v = D[f] + g.w[minmax(f, u)]
                v < best && ((best, arg) = (v, f))
            end
            arg > 0 && push!(new, (u, best, arg))
        end
        # the target cannot be reached: the target alone
        isempty(new) && return [target]
        foreach(((u, d, a),) -> (D[u] = d; pred[u] = a), new)
    end
    path = [target]
    while path[end] != start
        push!(path, pred[path[end]])
    end
    return reverse!(path)
end

# Lookup tables for _river_columns: a NamedTuple of two Dicts, `edge` (node
# pair as stored -> edge) and `diag` (quad diagonal as minmax -> column)
function _river_lookup(s)
    edge = Dict((s.edges[1, k], s.edges[2, k]) => k for k in axes(s.edges, 2))
    diag = Dict{Tuple{Int, Int}, Int}()
    for c in axes(s.corners, 2)
        ws, es, wn, en = s.corners[:, c]
        diag[minmax(ws, en)] = c; diag[minmax(es, wn)] = c
    end
    return (edge = edge, diag = diag)
end

# A node river as columns (co2lab's project_rivers_to_cells): along an edge
# the shallower column beside it (the only one at the boundary), across a
# quad diagonal the column it crosses; `look` from _river_lookup(s)
function _river_columns(s, look, nodes)
    edge, diag = look.edge, look.diag
    # repeated nodes dropped
    nodes = nodes[[true; diff(nodes) .!= 0]]
    cells = Int[]
    for (a, b) in zip(nodes[1:end-1], nodes[2:end])
        k = get(edge, minmax(a, b), 0)
        if k == 0
            push!(cells, diag[minmax(a, b)])
        else
            cols = sort(s.edge_columns[:, k])
            push!(cells, cols[1] == 0 ? cols[2] : cols[argmin(s.zc[cols])])
        end
    end
    return cells
end

# The part of a river outside traps, with one column into the traps at each
# end (co2lab's truncate_rivers_inside_traps; empty if it never leaves them)
function _truncate_river(river, traps)
    outside = findall(c -> traps[c] == 0, river)
    isempty(outside) && return Int[]
    return river[max(1, outside[1] - 1):min(length(river), outside[end] + 1)]
end
