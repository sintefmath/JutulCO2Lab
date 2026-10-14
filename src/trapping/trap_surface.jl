# The top surface as the trapping analysis sees it: what co2lab-spillpoint
# takes from the top-surface grid Gt -- node coordinates and depths, edges with
# the columns on either side, and per column its caprock depth, height, area,
# centroid and corner nodes.
#
# Nodes are numbered as co2lab numbers them (sorted by x, then y, then depth),
# so that ties between equally steep neighbours are broken the same way.
#
# Faults (columns' top faces not sharing an edge) are sealing by default: their
# edges are marked `sealed` and treated as closed boundaries. With open faults
# the surface is stitched as co2lab's topSurfaceGrid does (stitch_neighbors).

# ---------------------------------- Surface -----------------------------------

"""
    trap_surface(mesh, columns::VEColumns; faults = :sealing, sheet = :top)
        -> NamedTuple

The caprock of a VE domain, prepared for trap analysis: the polygons of
[`ve_top_surface`](@ref) with co2lab's node numbering, plus the quantities
co2lab-spillpoint reads from its top-surface grid `Gt`. This is the input of
[`trap_analysis`](@ref).

# Arguments

- `mesh`, `columns`: the 3D mesh and its column partition
  ([`ve_columns`](@ref)).
- `faults`: how faults are treated (see below); `:sealing` (default) or
  `:open`.
- `sheet`: which columns form the surface:
  - `:top` (default): the first VE block of every pillar, i.e. the
    formation's top; with one block per pillar, every column (as co2lab's
    `Gt`);
  - an integer: a sheet of [`ve_sheets`](@ref), a contiguous surface of
    connected top faces (e.g. under an internal layer). Its other boundary
    edges are open: CO2 reaching them leaves upwards, to the sheet above
    ([`sheet_spill_targets`](@ref)).

# Returns

A new NamedTuple, a superset of [`ve_top_surface`](@ref)'s result (so map
plots and [`ve_contour_lines`](@ref) accept it). Its "columns" are its
polygons, numbered `1:length(blocks)`; every per-column quantity here and in
[`trap_analysis`](@ref)'s result is per polygon.

- `x`, `y`, `z`: node coordinates and caprock depths (m), nodes sorted by
  `(x, y, z)` as co2lab's `uniqueNodes` (stitched nodes moved).
- `polygons[col]`: the column's top-face nodes.
- `corners`: `4 × ncol`, each column's west-south, east-south, west-north and
  east-north node (co2lab's `cellNodes`); 0 for columns without four distinct
  top-face nodes.
- `edges`: `2 × nedge` node pairs (smaller first) in co2lab's order. When
  every column has four corners this is `cartesian_face_reordering`'s: each
  column's west and east edges column by column, then all south, then all
  north edges; otherwise sorted.
- `edge_columns`: `2 × nedge`, the columns on either side (0 at the boundary).
- `sealed`: per edge, whether it is a boundary edge along a sealing fault.
- `zc`, `H`: caprock depth and height (m) per column
  ([`ve_column_geometry`](@ref)).
- `area`, `cx`, `cy`: area (m²) and centroid of each polygon in map view
  (co2lab's `Gt.cells.volumes`, `Gt.cells.centroids`).
- `xc`, `yc`: the mean of each polygon's nodes, as [`ve_top_surface`](@ref).
- `outline`: boundary edges as node pairs, as [`ve_top_surface`](@ref).
- `blocks`: the column of `columns` for each polygon (`1:n` with one block
  per pillar).

# Faults

A fault is a column side whose logical neighbour pillar holds a column of the
same surface, but which does not share its top edge with it. A fault without
throw leaves no gap in the grid and is not seen.

- `faults = :sealing`: closed boundaries. Their edges are marked in
  `sealed`; [`trap_analysis`](@ref) lets no CO2 spill across them.
- `faults = :open`: stitched as co2lab's `topSurfaceGrid` does, giving
  co2lab's surface: for columns connected in the 3D grid without a shared
  edge, the two closest node pairs (in map view) are merged at their mean
  position. Columns not connected in 3D (a complete offset) keep an open
  boundary between them, as in co2lab.

# Notes

- Each polygon has only the corners of its top face (at most four): nodes
  that faults add on the edges are left out, as co2lab's faces have none.
- Nodes at the same position are one node.
"""
function trap_surface(mesh, columns::VEColumns; faults::Symbol = :sealing,
                      sheet = :top)
    faults in (:sealing, :open) ||
        throw(ArgumentError("faults must be :sealing or :open, got :$faults"))
    blocks = _surface_blocks(mesh, columns, sheet)
    s = ve_top_surface(mesh, columns; blocks)
    geo = ve_column_geometry(mesh, columns)

    # corners only (faults leave extra nodes on edges), and the edges on faults
    keep = [_corner_positions(s.x, s.y, s.z, p) for p in s.polygons]
    at_fault = _fault_edges(mesh, columns, s, keep)
    polygons = [p[k] for (p, k) in zip(s.polygons, keep)]
    # co2lab's node order: sorted by x, y, z, equal positions one node
    x, y, z, polygons, renum = _unique_nodes(s.x, s.y, s.z, polygons)
    fault_pairs = Set(minmax(renum[a], renum[b]) for (a, b) in at_fault)
    # open faults: joined as co2lab's topSurfaceGrid does; nothing left to seal
    if faults == :open
        x, y, z, polygons = _stitch_faults(x, y, z, polygons,
                                           _column_connections(mesh, columns, blocks))
        empty!(fault_pairs)
    end

    edges, edge_columns = _surface_edges(polygons)
    corners = _polygon_corners(x, y, polygons)
    # co2lab's order on lattice-like surfaces (ties in the node-based method)
    if all(>(0), corners)
        order = _co2lab_edge_order(edges, corners)
        edges, edge_columns = edges[:, order], edge_columns[:, order]
    end
    # boundary edges along a fault
    sealed = [edge_columns[2, k] == 0 && Tuple(edges[:, k]) in fault_pairs
              for k in axes(edges, 2)]
    area, cx, cy = _polygon_geometry(x, y, polygons)
    xc = [sum(x[p]) / length(p) for p in polygons]
    yc = [sum(y[p]) / length(p) for p in polygons]
    return (x = x, y = y, z = z, polygons = polygons,
            corners = corners, edges = edges, edge_columns = edge_columns,
            sealed = sealed, zc = geo.z_top[blocks], H = geo.H[blocks],
            area = area, cx = cx, cy = cy, xc = xc, yc = yc,
            outline = _outline_edges(polygons), blocks = blocks)
end

# The nodes used by the polygons, sorted by (x, y, z), equal positions merged
# (co2lab's uniqueNodes): new coordinate vectors and polygons, and each old
# node's new number (0 if unused)
function _unique_nodes(x, y, z, polygons)
    pts = collect(zip(x, y, z))
    upts = sort!(unique(pts[unique(reduce(vcat, polygons))]))
    index = Dict(p => k for (k, p) in enumerate(upts))
    renum = [get(index, p, 0) for p in pts]
    return (first.(upts), getindex.(upts, 2), last.(upts),
            [renum[p] for p in polygons], renum)
end

# The positions of a top face's corners in its node list. A corner-point top
# face has at most four corners, but faults beside it add nodes on its edges
# (co2lab's faces have none); the node nearest the line through its two
# neighbours is dropped until four distinct nodes are left.
function _corner_positions(x, y, z, p)
    pos = collect(eachindex(p))
    p = copy(p)
    while length(unique(p)) > 4
        dev = map(eachindex(p)) do k
            a, b, c = p[mod1(k - 1, length(p))], p[k], p[mod1(k + 1, length(p))]
            u = (x[c] - x[a], y[c] - y[a], z[c] - z[a])
            v = (x[b] - x[a], y[b] - y[a], z[b] - z[a])
            # distance of b from the line a-c, relative to |a - c|
            w = (u[2]v[3] - u[3]v[2], u[3]v[1] - u[1]v[3], u[1]v[2] - u[2]v[1])
            sqrt(sum(abs2, w)) / sum(abs2, u)
        end
        k = argmin(dev)
        deleteat!(p, k); deleteat!(pos, k)
    end
    return pos
end

# -------------------------------- Fault edges ---------------------------------

# The edges of the corner polygons that lie on a fault, as a set of node pairs
# of `s` (ve_top_surface). An edge is on a fault if it lies on a column side
# whose logical neighbour pillar holds a column of the surface; the side is
# found from the top cell's lateral faces. `keep`: each polygon's corner
# positions.
function _fault_edges(mesh, columns, s, keep)
    ijk = [cell_ijk(mesh, c) for c in 1:number_of_cells(mesh)]
    bdir = _boundary_directions(mesh)
    pillars = Set(columns.column_ij[s.blocks])
    out = Set{Tuple{Int, Int}}()
    for (i, p) in enumerate(s.polygons)
        col = s.blocks[i]
        top = columns.cells[columns.pos[col]]
        sides = _lateral_face_nodes(mesh, ijk, bdir, top)
        fault(d) = _logical_neighbour_pillar(columns.column_ij[col], d) in pillars
        k, n = keep[i], length(p)
        for (i, a) in enumerate(k)
            b = k[mod1(i + 1, length(k))]
            # the raw pieces from corner a to corner b, cyclically
            chain = [p[mod1(t, n)] for t in a:(b > a ? b : b + n)]
            on_fault = any(zip(chain[1:end-1], chain[2:end])) do (u, v)
                d = _edge_side(sides, s.mesh_nodes[u], s.mesh_nodes[v])
                d > 0 && fault(d)
            end
            on_fault && push!(out, minmax(p[a], p[b]))
        end
    end
    return out
end

# The lateral faces of cell c as (direction, nodes) — boundary faces by their
# tag, interior faces by the logical position of the cell beyond
function _lateral_face_nodes(mesh, ijk, bdir, c)
    out = Tuple{Int, Vector{Int}}[]
    for f in mesh.boundary_faces.cells_to_faces[c]
        nodes = collect(mesh.boundary_faces.faces_to_nodes[f])
        bdir[f] <= 4 && push!(out, (bdir[f], nodes))
    end
    for f in mesh.faces.cells_to_faces[c]
        d, _ = _interior_face_direction(mesh, ijk, c, f)
        d <= 4 && push!(out, (d, collect(mesh.faces.faces_to_nodes[f])))
    end
    return out
end

# The direction of the lateral face holding mesh nodes u and v (0: none)
function _edge_side(sides, u, v)
    u == v && return 0
    i = findfirst(((d, nodes),) -> u in nodes && v in nodes, sides)
    return isnothing(i) ? 0 : first(sides[i])
end

# ------------------------------- Fault stitching ------------------------------

# Pairs (a < b) of the surface's polygons whose columns (`blocks`) are
# connected by a face of the 3D grid
function _column_connections(mesh, columns, blocks)
    # each cell's polygon (0: its column is not on the surface)
    poly = zeros(Int, number_of_columns(columns))
    poly[blocks] = eachindex(blocks)
    col = [b == 0 ? 0 : poly[b] for b in columns.partition]
    pairs = Set{Tuple{Int, Int}}()
    for (c1, c2) in eachcol(get_neighborship(mesh))
        a, b = col[c1], col[c2]
        a > 0 && b > 0 && a != b && push!(pairs, minmax(a, b))
    end
    return pairs
end

# co2lab's stitch_surface_discontinuities: for each pair of connected columns
# without a shared edge, the two closest node pairs in map view (first of
# equal distances, in co2lab's order) are merged — each group of merged nodes
# at its mean position — and the nodes renumbered
function _stitch_faults(x, y, z, polygons, connections)
    # boundary edges have columns (c, 0)
    shared = Set(minmax(c[1], c[2]) for c in eachcol(_surface_edges(polygons)[2])
                 if c[2] > 0)
    merge = Tuple{Int, Int}[]
    for (a, b) in sort!(collect(setdiff(connections, shared)))
        n1, n2 = _edge_node_order(polygons[a]), _edge_node_order(polygons[b])
        combos = [(i, j) for i in n1 for j in n2]
        d = [(x[i] - x[j])^2 + (y[i] - y[j])^2 for (i, j) in combos]
        for k in sortperm(d)[1:2]
            i, j = combos[k]
            i != j && push!(merge, minmax(i, j))
        end
    end
    isempty(merge) && return (x, y, z, polygons)

    x, y, z = copy(x), copy(y), copy(z)
    for group in _node_groups(length(x), merge)
        w = 1 / length(group)
        # co2lab: a sparse row of 1/n times the coordinates, in node order
        x[group] .= foldl((s, n) -> s + w * x[n], group; init = 0.0)
        y[group] .= foldl((s, n) -> s + w * y[n], group; init = 0.0)
        z[group] .= foldl((s, n) -> s + w * z[n], group; init = 0.0)
    end
    return _unique_nodes(x, y, z, polygons)[1:4]
end

# A polygon's nodes as co2lab lists them (padded_nodelist): each edge's nodes,
# smaller index first, in edge order, repeats removed
_edge_node_order(p) =
    unique([n for (k, a) in enumerate(p) for n in minmax(a, p[mod1(k + 1, length(p))])])

# Groups (sorted node lists) of nodes joined by the pairs `merge`
function _node_groups(nn, merge)
    g = Graphs.SimpleGraph(nn)
    foreach(((a, b),) -> Graphs.add_edge!(g, a, b), merge)
    return [sort(c) for c in Graphs.connected_components(g) if length(c) > 1]
end

# ------------------------------- Edges, corners -------------------------------

# Edges as sorted node pairs (2 × ne) and the columns on either side (2 × ne,
# the smaller column first; 0 for a boundary edge)
function _surface_edges(polygons)
    # the columns on each edge (one at the boundary); collapsed edges skipped
    sides = Dict{Tuple{Int, Int}, Vector{Int}}()
    for (col, p) in enumerate(polygons), (k, a) in enumerate(p)
        b = p[k == length(p) ? 1 : k + 1]
        a == b && continue
        push!(get!(sides, minmax(a, b), Int[]), col)
    end
    pairs = sort!(collect(keys(sides)))
    edges = [e[i] for i in 1:2, e in pairs]
    cols = [i <= length(sides[e]) ? sort(sides[e], rev = true)[i] : 0
            for i in 1:2, e in pairs]
    # boundary edges: (column, 0); interior: (smaller, larger)
    for k in axes(cols, 2)
        cols[2, k] != 0 && (cols[:, k] = sort(cols[:, k]))
    end
    return (edges, cols)
end

# co2lab's edge order on a lattice (topSurfaceGrid's cartesian_face_reordering):
# west and east edge of each column in turn, then the south edges, then the
# north edges, each edge where it first appears
function _co2lab_edge_order(edges, corners)
    index = Dict((edges[1, k], edges[2, k]) => k for k in axes(edges, 2))
    e(a, b) = index[minmax(a, b)]
    ws, es, wn, en = (corners[i, :] for i in 1:4)
    # each column's west and east edges in turn, then all south, then all north
    first = [k for c in axes(corners, 2) for k in (e(ws[c], wn[c]), e(es[c], en[c]))]
    second = [[e(ws[c], es[c]) for c in axes(corners, 2)];
              [e(wn[c], en[c]) for c in axes(corners, 2)]]
    order = unique([unique(first); unique(second)])
    # edges left over: not a lattice, keep the sorted order
    length(order) == size(edges, 2) || return collect(axes(edges, 2))
    return order
end

# West-south, east-south, west-north, east-north node of each column
# (4 × ncol; 0 where the polygon does not have four distinct nodes)
function _polygon_corners(x, y, polygons)
    corners = zeros(Int, 4, length(polygons))
    for (col, p) in enumerate(polygons)
        q = unique(p)
        length(q) == 4 || continue
        # the two southern nodes, then the two northern, each west to east
        s = sort(q, by = n -> (y[n], x[n]))
        south, north = sort(s[1:2], by = n -> x[n]), sort(s[3:4], by = n -> x[n])
        corners[:, col] = [south; north]
    end
    return corners
end

# Area and centroid of each polygon in map view (shoelace formula, relative to
# the polygon's first node: map coordinates of ~1e6 m would cancel badly)
function _polygon_geometry(x, y, polygons)
    n = length(polygons)
    area, cx, cy = zeros(n), zeros(n), zeros(n)
    for (col, p) in enumerate(polygons)
        x0, y0 = x[p[1]], y[p[1]]
        a, sx, sy = 0.0, 0.0, 0.0
        for (k, i) in enumerate(p)
            j = p[k == length(p) ? 1 : k + 1]
            xi, yi, xj, yj = x[i] - x0, y[i] - y0, x[j] - x0, y[j] - y0
            w = xi * yj - xj * yi
            a += w
            sx += (xi + xj) * w
            sy += (yi + yj) * w
        end
        area[col], cx[col], cy[col] = abs(a) / 2, x0 + sx / (3a), y0 + sy / (3a)
    end
    return (area, cx, cy)
end
