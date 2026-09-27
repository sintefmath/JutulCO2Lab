# The top surface as the trapping analysis sees it (roadmap Stage 8 item 1) —
# what co2lab-spillpoint takes from the top-surface grid Gt: node coordinates
# and depths, edges with the columns on either side, and per column its
# caprock depth, height, area, centroid and corner nodes.
#
# Nodes are numbered as co2lab numbers them (sorted by x, then y, then depth),
# so that ties between equally steep neighbours are broken the same way.
#
# Across faults the columns' top faces do not share edges; co2lab's
# topSurfaceGrid stitches the surface there (stitch_neighbors): for each pair
# of columns connected in the 3D grid but not along an edge, the two closest
# node pairs in map view are merged into their mean. The same is done here.

# ---------------------------------- Surface -----------------------------------

"""
    trap_surface(mesh, columns::VEColumns) -> NamedTuple

The caprock of a VE domain for trap analysis: [`ve_top_surface`](@ref)'s
polygons with co2lab's numbering and the quantities co2lab-spillpoint reads
from its top-surface grid `Gt`.

Each column's polygon has only the corners of its top face (at most four;
nodes that faults add on the edges are left out, as co2lab's faces have
none). The surface is continuous: nodes at the same position are one node, and
across faults (columns connected in the 3D grid without a shared edge) it
is stitched as co2lab's `topSurfaceGrid` stitches it — the two closest node
pairs of the two columns, in map view, merged at their mean position.

Returns the NamedTuple `(x, y, z, polygons, corners, edges, edge_columns, zc,
H, area, cx, cy, xc, yc, outline)`, all new — a superset of
[`ve_top_surface`](@ref)'s result, so that the map plots and
[`ve_contour_lines`](@ref) take it too:
- `x`, `y`, `z`: node coordinates and caprock depths, nodes sorted by
  `(x, y, z)` as co2lab's `uniqueNodes` (stitched nodes moved);
- `polygons[col]`: the column's top-face nodes;
- `corners`: `4 × ncol`, each column's west-south, east-south, west-north and
  east-north nodes (co2lab's `cellNodes`; for columns with four distinct
  top-face nodes, else 0);
- `edges`: `2 × nedge` node pairs (smaller first) in co2lab's order — when
  every column has four corners, as its `cartesian_face_reordering`: each
  column's west and east edges column by column, then all south edges,
  then all north edges; otherwise sorted — and
  `edge_columns`: `2 × nedge`, the columns on either side (0 at the
  boundary);
- `zc`, `H`: caprock depth and height per column
  ([`ve_column_geometry`](@ref));
- `area`, `cx`, `cy`: the area and centroid of each column's polygon in map
  view (co2lab's `Gt.cells.volumes`, `Gt.cells.centroids`);
- `xc`, `yc`: the mean of each polygon's nodes, as [`ve_top_surface`](@ref);
- `outline`: boundary edges as node pairs, as [`ve_top_surface`](@ref).
`mesh` and `columns` are only read.
"""
function trap_surface(mesh, columns::VEColumns)
    s = ve_top_surface(mesh, columns)
    geo = ve_column_geometry(mesh, columns)

    # corners only (faults leave extra nodes on edges); co2lab's node order:
    # sorted by x, y, z, equal positions one node
    polygons = [_corner_nodes(s.x, s.y, s.z, p) for p in s.polygons]
    x, y, z, polygons = _unique_nodes(s.x, s.y, s.z, polygons)
    x, y, z, polygons = _stitch_faults(x, y, z, polygons,
                                       _column_connections(mesh, columns))

    edges, edge_columns = _surface_edges(polygons)
    corners = _polygon_corners(x, y, polygons)
    # co2lab's order on lattice-like surfaces (ties in the node-based method)
    if all(>(0), corners)
        order = _co2lab_edge_order(edges, corners)
        edges, edge_columns = edges[:, order], edge_columns[:, order]
    end
    area, cx, cy = _polygon_geometry(x, y, polygons)
    xc = [sum(x[p]) / length(p) for p in polygons]
    yc = [sum(y[p]) / length(p) for p in polygons]
    return (x = x, y = y, z = z, polygons = polygons,
            corners = corners, edges = edges,
            edge_columns = edge_columns, zc = geo.z_top, H = geo.H,
            area = area, cx = cx, cy = cy, xc = xc, yc = yc,
            outline = _outline_edges(polygons))
end

# The nodes used by the polygons, sorted by (x, y, z), equal positions merged
# (co2lab's uniqueNodes): new coordinate vectors and polygons
function _unique_nodes(x, y, z, polygons)
    pts = collect(zip(x, y, z))
    upts = sort!(unique(pts[unique(reduce(vcat, polygons))]))
    index = Dict(p => k for (k, p) in enumerate(upts))
    return (first.(upts), getindex.(upts, 2), last.(upts),
            [[index[pts[n]] for n in p] for p in polygons])
end

# A top face's corner nodes: a corner-point cell's top face has at most four;
# where a fault splits the faces beside it, GeoEnergyIO adds nodes on its
# edges. The node closest to the line through its two neighbours is dropped,
# repeatedly, until four distinct nodes are left (co2lab's faces have none).
function _corner_nodes(x, y, z, p)
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
        deleteat!(p, argmin(dev))
    end
    return p
end

# ------------------------------- Fault stitching ------------------------------

# Pairs (a < b) of columns connected by a face of the 3D grid
function _column_connections(mesh, columns)
    col = columns.partition
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
    return _unique_nodes(x, y, z, polygons)
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
    first = [k for c in axes(corners, 2) for k in (e(ws[c], wn[c]), e(es[c], en[c]))]
    second = [[e(ws[c], es[c]) for c in axes(corners, 2)];
              [e(wn[c], en[c]) for c in axes(corners, 2)]]
    order = unique([unique(first); unique(second)])
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
