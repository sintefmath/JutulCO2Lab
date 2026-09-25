# VE reservoir domain on a CoarseMesh, and its lateral boundary sides.
#
# See roadmap Stage 1 and analysis/ve-meshing/grid-representation-comparison.md.
# Everything is derived from the 3D mesh and the VEColumns partition.

"""
    ve_average_rock(columns, dz; perm_x, poro, ntg = nothing) -> NamedTuple

Thickness-weighted column averages, exactly as MRST's `averageRock`:
`sum(x .* dz) / sum(dz)` per column, for x-permeability, porosity and (if
given) net-to-gross. `dz` is per fine cell in `columns.cells` order (from
[`ve_column_geometry`](@ref)); `perm_x`, `poro`, `ntg` are per fine cell of the
3D mesh. Only x-permeability is used — co2lab discards anisotropy.

Returns the NamedTuple `(perm, poro, ntg)` of new per-column vectors;
`ntg` is `nothing` if not given. The arguments are only read.
"""
function ve_average_rock(columns::VEColumns, dz::AbstractVector;
                         perm_x, poro, ntg = nothing)
    avg(x) = [sum(x[columns.cells[r]] .* dz[r]) / sum(dz[r])
              for r in (columns.pos[c]:(columns.pos[c+1] - 1)
                        for c in 1:number_of_columns(columns))]
    return (perm = avg(perm_x), poro = avg(poro),
            ntg = isnothing(ntg) ? nothing : avg(ntg))
end

"""
    ve_domain(mesh, columns; perm_x, poro, ntg = nothing) -> NamedTuple

Build the VE reservoir domain: a `CoarseMesh` over `mesh`, one block per
column, with thickness-averaged rock and the gravity term taken from the
caprock.

Returns the NamedTuple `(domain, coarse, geometry, rock)`, all new:
- `domain`: the `DataDomain` on the `CoarseMesh` (from `reservoir_domain`),
  with the averaged `:permeability`, `:porosity` (and `:net_to_gross`) per
  column and `:gdz` on `Faces()` from the caprock depths;
- `coarse`: that `CoarseMesh` — the same object `domain` wraps;
- `geometry`: [`ve_column_geometry`](@ref)'s result;
- `rock`: [`ve_average_rock`](@ref)'s result.
`mesh` and `columns` are only read.

Deliberately **not** done, and why:

- Permeability is **not** multiplied by the column height. On a `CoarseMesh`
  the height is already in the geometry (face area ≈ edge length × H), so
  Jutul's default transmissibility reproduces MRST's `perm·H` one; scaling
  permeability too would double-count.
- Cell centroids are **left mid-column**. Transmissibilities are computed from
  them; moving their depth to the caprock would change every transmissibility.
  The caprock datum enters through `:gdz` instead, which takes precedence over
  the centroid-based default (`JutulDarcy/src/multiphase.jl:303-308`). Anything
  else that reads centroid depths must be overridden where it is used — see
  roadmap Stage 1 item 7.
"""
function ve_domain(mesh, columns::VEColumns; perm_x, poro, ntg = nothing)
    geometry = ve_column_geometry(mesh, columns)
    rock = ve_average_rock(columns, geometry.dz; perm_x, poro, ntg)
    coarse = CoarseMesh(mesh, columns.partition)
    kw = isnothing(rock.ntg) ? (;) : (net_to_gross = rock.ntg,)
    domain = reservoir_domain(coarse; permeability = rock.perm, porosity = rock.poro,
                              kw...)
    domain[:gdz, Faces()] = Jutul.compute_face_gdz(domain[:neighbors], geometry.z_top)
    return (domain = domain, coarse = coarse, geometry = geometry, rock = rock)
end

"""
    VEBoundarySides

Lateral boundary sides of the VE columns — the counterpart of the boundary
edges of MRST's top-surface grid. One entry per (column, side):

- `column`: the column (coarse cell)
- `edge_depth`: caprock depth along the side — mean depth of the two top-face
  nodes bounding it (MRST's `Gt.faces.z`, but without MRST's averaging of
  nodes across pillar offsets)
- `edge_midpoint`: (x, y) midpoint of that top edge
- `area`, `centroid`, `normal`: aggregate of the column's lateral boundary
  fine faces on this side (outward unit normal)
- `half_trans`: two-point half transmissibility from the column's centroid to
  the side, `k·A·(c·n)/|c|²` — the same formula family as the `CoarseMesh`'s
  interior transmissibilities
- `fine_faces`: boundary-face indices in the 3D mesh
- `across_fault`: a column exists at the logical neighbour position, but no
  fine face connects to it — a fault with complete offset (Johansen: 88 sides,
  along its main fault and one more). Such a side is an outer side, as in
  MRST's top-surface grid, but physically sealed: give it a boundary
  condition only deliberately (MRST's Johansen example keeps its main fault
  sealing, `makeJohansenVEgrid.m:143-147`).
"""
struct VEBoundarySides
    column::Vector{Int}
    edge_depth::Vector{Float64}
    edge_midpoint::Vector{NTuple{2, Float64}}
    area::Vector{Float64}
    centroid::Vector{NTuple{3, Float64}}
    normal::Vector{NTuple{3, Float64}}
    half_trans::Vector{Float64}
    fine_faces::Vector{Vector{Int}}
    across_fault::Vector{Bool}
end

Base.length(s::VEBoundarySides) = length(s.column)

"""
    ve_boundary_sides(mesh, columns, domain) -> (sides, unassigned)

The outer boundary sides of the columns: the column sides
([`ve_column_sides`](@ref)) with no neighbouring column. Everything is
decided by topology — face directions and cell connectivity — and geometry
only supplies the quantities.

Where a column has a neighbour on some side but stands higher (an internal
caprock step), the exposed part of that side consists of boundary faces of
the 3D grid, but the side is not an outer boundary: MRST's top-surface grid
has an interior edge there (on Johansen: 173 such cases; the remaining 456
sides are exactly MRST's 456 boundary edges). Such faces are returned in
`unassigned`; they stay no-flow, as in MRST.

The edge of a side is where the column's top face meets the top cell's
boundary faces on that side: exactly two shared nodes. `mesh`, `columns` and
`domain` are only read.

Returns `(sides, unassigned)`:
- `sides`: a new [`VEBoundarySides`](@ref), one entry per outer side, ordered
  by column and, within a column, by direction (left, right, lower, upper);
- `unassigned`: `Vector{Int}` of boundary-face indices of `mesh` on sides that
  do have a neighbour (exposed fragments; no-flow), in the same order.
"""
function ve_boundary_sides(mesh::UnstructuredMesh, columns::VEColumns, domain)
    cs = ve_column_sides(mesh, columns)
    g = _boundary_geometry(mesh)
    bdir = _boundary_directions(mesh)
    perm, cc = domain[:permeability], domain[:cell_centroids]
    out = (column = Int[], depth = Float64[], mid = NTuple{2, Float64}[],
           area = Float64[], cen = NTuple{3, Float64}[], nrm = NTuple{3, Float64}[],
           T = Float64[], faces = Vector{Int}[], fault = Bool[])
    column_at = Dict(columns.pillar[c] => c for c in 1:number_of_columns(columns))
    unassigned = Int[]
    for s in 1:length(cs)
        col, fs = cs.column[s], cs.boundary_faces[s]
        if cs.neighbour[s] != 0
            append!(unassigned, fs); continue
        end
        d, top = cs.direction[s], first(column_cells(columns, col))
        depth, mid = _side_top_edge(g, bdir, top, d, col)
        centroid = (cc[1, col], cc[2, col], cc[3, col])
        A, C, N, t = _side_aggregate(g, fs, perm[col], centroid)
        push!(out.column, col); push!(out.depth, depth); push!(out.mid, mid)
        push!(out.area, A); push!(out.cen, C); push!(out.nrm, N)
        push!(out.T, t); push!(out.faces, fs)
        # a column exists across this side, but nothing connects to it: sealed fault
        across = _logical_neighbour_pillar(columns.pillar[col], d)
        push!(out.fault, haskey(column_at, across))
    end
    sides = VEBoundarySides(out.column, out.depth, out.mid, out.area, out.cen, out.nrm,
                            out.T, out.faces, out.fault)
    return (sides, unassigned)
end

# Logical (i, j) of the neighbouring pillar in lateral direction d (1..4 =
# left, right, lower, upper).
_logical_neighbour_pillar(ij, d) = ij .+ ((-1, 0), (1, 0), (0, -1), (0, 1))[d]

# Boundary-face geometry of the 3D mesh, bundled for the helpers below.
function _boundary_geometry(mesh::UnstructuredMesh)
    geo = tpfv_geometry(mesh)
    bf = mesh.boundary_faces
    return (bn = geo.boundary_normals, bc = geo.boundary_centroids,
            ba = geo.boundary_areas, c2f = bf.cells_to_faces, nodes = bf.faces_to_nodes,
            X = mesh.node_points)
end

"""
    _side_top_edge(g, bdir, top, d, col) -> (depth, midpoint)

Caprock edge of column `col`'s side in direction `d`: the two nodes shared by
the top cell's top face and its boundary faces on that side. Returns their
mean depth and (x, y) midpoint.
"""
function _side_top_edge(g, bdir, top, d, col)
    tops = [f for f in g.c2f[top] if bdir[f] == _TOP]
    length(tops) == 1 ||
        error("column $col: top cell $top has $(length(tops)) top boundary faces, " *
              "expected 1")
    side_nodes = Set{Int}()
    for f in g.c2f[top]
        bdir[f] == d && union!(side_nodes, g.nodes[f])
    end
    shared = [n for n in g.nodes[tops[1]] if n in side_nodes]
    length(shared) == 2 ||
        error("column $col, $(_DIRECTIONS[d]) side: top edge has " *
              "$(length(shared)) nodes, expected 2")
    Xa, Xb = g.X[shared[1]], g.X[shared[2]]
    return ((Xa[3] + Xb[3]) / 2, ((Xa[1] + Xb[1]) / 2, (Xa[2] + Xb[2]) / 2))
end

# Aggregate one side's fine faces fs: area, centroid, outward unit normal, and
# the two-point half transmissibility k·A·(c·n)/|c|² from the column centroid.
function _side_aggregate(g, fs, perm, col_centroid)
    A = sum(g.ba[f] for f in fs)
    C = ntuple(d -> sum(g.ba[f] * g.bc[d, f] for f in fs) / A, 3)
    Nv = ntuple(d -> sum(g.ba[f] * g.bn[d, f] for f in fs), 3)
    N = Nv ./ sqrt(sum(abs2, Nv))
    cvec = C .- col_centroid
    t = perm * A * sum(cvec .* N) / sum(abs2, cvec)
    return (A, C, N, t)
end
