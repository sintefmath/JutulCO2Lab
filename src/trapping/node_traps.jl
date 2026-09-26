# Structural traps along the surface's nodes and edges: co2lab-spillpoint's
# node-based method, trapAnalysis(Gt, false) — findUpslopeNeighbor,
# nodeSpillField, computeNodeTraps (roadmap Stage 9 item 3). Results per node;
# node_to_cells.jl projects them to columns (n2cTraps).
#
# Every node drains to its steepest upslope neighbour (over edges and quad
# diagonals); nodes draining to the same local top form a region; a region
# spills over the shallowest point of its rim; regions whose spill paths form
# a cycle are merged into one trap, cycle by cycle, as co2lab does.

# ------------------------------ Upslope neighbour -----------------------------

# Node neighbours over edges (sorted), leaving out edges touching `closed`
function _node_neighbours(s, closed)
    nn = length(s.z)
    nbrs = [Int[] for _ in 1:nn]
    for k in axes(s.edges, 2)
        a, b = s.edges[1, k], s.edges[2, k]
        (a in closed || b in closed) && continue
        push!(nbrs[a], b); push!(nbrs[b], a)
    end
    return map(v -> sort!(unique(v)), nbrs)
end

# co2lab's findUpslopeNeighbor: among a node, its edge neighbours and the
# nodes sharing exactly two neighbours with it (quad diagonals), the one with
# the largest signed squared slope sign(dz) dz² / |d|² (d the 3D offset);
# the first of equal maxima in node order
function _upslope_neighbours(s, closed)
    nbrs = _node_neighbours(s, closed)
    up = zeros(Int, length(s.z))
    for n in eachindex(s.z)
        common = Dict{Int, Int}()
        for k in nbrs[n], j in nbrs[k]
            j == n || (common[j] = get(common, j, 0) + 1)
        end
        cand = sort!(unique([n; nbrs[n]; [j for (j, c) in common if c == 2]]))
        best = -Inf
        for c in cand
            dx, dy, dz = s.x[n] - s.x[c], s.y[n] - s.y[c], s.z[n] - s.z[c]
            slope = sign(dz) * dz^2 / (dx^2 + dy^2 + dz^2)
            isnan(slope) && (slope = 0.0)
            slope > best && ((best, up[n]) = (slope, c))
        end
    end
    return up
end

# --------------------------------- Regions ------------------------------------

# co2lab's nodeSpillField: the region of each node (0: drains out of the
# domain over an open boundary node; k: to the k-th local top in node order),
# the downstream neighbour (0 at tops and open boundary nodes) and the spill
# edges — rows (edge, region 1, region 2, z1, z2, node 1, node 2) of the
# edges between regions, in edge order
function _spill_field(s, up, closed_bnodes)
    nn = length(s.z)
    boundary = fill(false, nn)
    for k in axes(s.edges, 2)
        any(==(0), s.edge_columns[:, k]) && (boundary[s.edges[:, k]] .= true)
    end
    open = boundary .& .!in.(1:nn, Ref(Set(closed_bnodes)))
    tops = [up[n] == n && !open[n] for n in 1:nn]

    region = fill(-1, nn)
    region[open] .= 0
    region[tops] .= 1:count(tops)
    for n in 1:nn
        # follow the upslope chain to a node with a region
        chain = Int[]
        m = n
        while region[m] < 0
            push!(chain, m); m = up[m]
            length(chain) > nn && error("upslope chain without end at node $n")
        end
        region[chain] .= region[m]
    end

    down = copy(up)
    down[open .| tops] .= 0
    rows = [(k, region[a], region[b], s.z[a], s.z[b], a, b)
            for (k, (a, b)) in enumerate(eachcol(s.edges)) if region[a] != region[b]]
    return (region = region, down = down, spill = rows)
end

# ------------------------- Spill point of a region group -----------------------

# co2lab's locate_neighbor_traps: for a group of regions `members`, its spill
# depth (the shallowest of its rim edges' deeper ends), the spill edges at
# that depth (one per inner node: the one with the shallowest outer end) and
# the groups beyond them (`group_of` maps region → group)
function _group_spill(members, group_of, spill)
    member_set = Set(members)
    inside(r) = r in member_set
    rim = [(i, max(e[4], e[5])) for (i, e) in enumerate(spill)
           if inside(e[2]) != inside(e[3])]
    # stable sort by depth, as co2lab's sortrows
    sort!(rim, by = last)
    z_spill = rim[1][2]
    tied = [i for (i, z) in rim if z == z_spill]

    # per inner node, the tied edge with the shallowest outer end: co2lab lists
    # the edges with the group on side 1 first, then on side 2, and sorts
    # stably by outer depth, then by node
    ori = [[(spill[i][6], spill[i][5], i) for i in tied if inside(spill[i][2])];
           [(spill[i][7], spill[i][4], i) for i in tied if inside(spill[i][3])]]
    ori = sort(sort(ori, by = r -> r[2]), by = r -> r[1])
    chosen = [ori[findfirst(r -> r[1] == n, ori)][3] for n in unique(first.(ori))]

    beyond = setdiff(unique([r for i in chosen for r in spill[i][2:3]]), [members; 0])
    groups = sort!(unique([group_of[r] for r in beyond]))
    return (groups = groups, z = z_spill, edges = [spill[i][1] for i in chosen])
end

# ------------------------------- Cycle merging --------------------------------

# 1 where A is nonzero (co2lab's spones)
_spones(A) = (B = copy(A); nonzeros(B) .= 1.0; dropzeros!(B))

# co2lab's merge_rows_cols: rows and columns of `cycle` ORed into the first,
# the others removed
function _merge_groups(A, cycle)
    n = size(A, 1)
    keep = setdiff(1:n, cycle[2:end])
    new = [c in cycle ? findfirst(==(cycle[1]), keep) : findfirst(==(c), keep)
           for c in 1:n]
    P = sparse(1:n, new, 1.0, n, length(keep))
    return _spones(P' * A * P)
end

# co2lab's track_stream: groups of regions (initially one each) and which
# group spills into which; walks of growing length are accumulated, and a
# group on a cycle is merged with every group on its cycles — one cycle at a
# time, the merged group spilling anew over its own shallowest rim point.
# Returns (members, z, spill_edges, adjacency), groups ordered by their first
# region.
function _merge_cycles(nreg, spill)
    members = [[r] for r in 1:nreg]
    group_of = collect(1:nreg)
    cmat = spzeros(nreg, nreg)
    z, edges = zeros(nreg), Vector{Vector{Int}}(undef, nreg)
    for g in 1:nreg
        r = _group_spill(members[g], group_of, spill)
        cmat[g, r.groups] .= 1.0
        z[g], edges[g] = r.z, r.edges
    end
    icmat, tmp = spzeros(nreg, nreg), copy(cmat)
    while nnz(tmp) > 0
        icmat = _spones(icmat + tmp)
        while (i = findfirst(!iszero, diag(icmat))) !== nothing
            cycle = [k for k in axes(icmat, 1) if icmat[i, k] != 0 && icmat[k, i] != 0]
            g = cycle[1]
            cmat, icmat, tmp = (_merge_groups(M, cycle) for M in (cmat, icmat, tmp))
            cmat[g, g] = 0.0; icmat[g, g] = 0.0
            dropzeros!(cmat); dropzeros!(icmat)

            members[g] = sort!(vcat(members[cycle]...))
            for v in (members, z, edges)
                deleteat!(v, cycle[2:end])
            end
            foreach(((k, m),) -> group_of[m] .= k, enumerate(members))
            r = _group_spill(members[g], group_of, spill)
            cmat[g, r.groups] .= 1.0
            z[g], edges[g] = r.z, r.edges
        end
        tmp = _spones(tmp * cmat)
    end
    return (members = members, z = z, spill_edges = edges, adjacency = cmat)
end

# ------------------------------- Node traps -----------------------------------

# co2lab's computeNodeTraps: trap and spill region of each node, the traps'
# spill depths and connections, and the rivers (node lists) from each trap's
# spill edges, following the downstream neighbours from both ends
function _node_traps(s, closed_bnodes, closed_fnodes)
    up = _upslope_neighbours(s, Set(closed_fnodes))
    f = _spill_field(s, up, closed_bnodes)
    g = _merge_cycles(maximum(f.region; init = 0), f.spill)

    nn, nt = length(s.z), length(g.members)
    trap_regions, traps = zeros(Int, nn), zeros(Int, nn)
    group_of = Dict(r => t for (t, m) in enumerate(g.members) for r in m)
    for n in 1:nn
        f.region[n] == 0 && continue
        t = group_of[f.region[n]]
        trap_regions[n] = t
        s.z[n] <= g.z[t] && (traps[n] = t)
    end
    rivers = [[_node_river(s, f.down, traps, t, e) for e in g.spill_edges[t]]
              for t in 1:nt]
    return (traps = traps, trap_z = g.z, trap_regions = trap_regions,
            trap_adj = g.adjacency .!= 0, rivers = rivers, down = f.down)
end

# The river over spill edge e of trap t: both ends extended downstream, then
# oriented to start in the trap (co2lab's trace_rivers)
function _node_river(s, down, traps, t, e)
    river = [s.edges[1, e], s.edges[2, e]]
    while down[river[end]] != 0
        push!(river, down[river[end]])
    end
    while down[river[1]] != 0
        pushfirst!(river, down[river[1]])
    end
    traps[river[1]] == t || traps[river[end]] == t ||
        error("river over edge $e does not start in trap $t")
    traps[river[1]] == t || reverse!(river)
    return river
end
