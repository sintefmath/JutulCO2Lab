# The faces of a hybrid VE grid: each coarse face's category, which decides
# how flux across it is computed, and the data that flux needs.
#
# Categories:
#   0 sealed: all its fine faces sealed (multiplier 0); no flow.
#   1 fine–fine: the ordinary 3D flux.
#   2 VE–VE, continuous: the two blocks' top cells share a face (the caprock
#     continues) and the columns face each other over their whole height;
#     the ordinary whole-column VE flux.
#   3 VE–VE, discontinuous: lateral, but the top cells are not connected (a
#     fault offset larger than the top cell, a pinch-out, stacked blocks), or
#     a block faces only part of the other (a short column above a layer
#     beside a full one). Virtual cell, each side over its interval at the
#     face.
#   4 VE–VE, vertical: stacked blocks through a leaky layer. Virtual cell over
#     the blocks' adjacent cells.
#   5 VE–fine, lateral or vertical: virtual cell over the fine cell's extent.
#
# Categories 3–5 are the transitions. They carry, per side, the depth interval
# of that side's cells at the face (below its block's top), a common depth
# (the fine cell's centroid for 5, the faces' mean depth for 3 and 4), and the
# transmissibility summed from the fine grid's, each fine face's times its
# multiplier. Category 2 keeps the coarse mesh's transmissibility (as pure VE),
# times its fine faces' area-weighted mean multiplier; 1 takes the fine grid's
# (with its multiplier).
#
# `fault_flux = :virtual` (default) follows the topology. `:co2lab` makes
# lateral VE–VE faces between pillars with one VE block each continuous
# (co2lab's smooth flux across faults, for comparisons with co2lab); stacked
# blocks follow the topology in any case.

const _FACE_SEALED = 0                 # all fine faces sealed: no flow
const _FACE_FINE = 1                   # fine–fine: 3D flux
const _FACE_VE = 2                     # continuous VE–VE: whole-column VE flux
const _FACE_VE_DISCONTINUOUS = 3       # lateral VE–VE, caprock not continuous
const _FACE_VE_VERTICAL = 4            # stacked VE–VE through a leaky layer
const _FACE_VE_FINE = 5                # VE–fine
const _TRANSITIONS = (_FACE_VE_DISCONTINUOUS, _FACE_VE_VERTICAL, _FACE_VE_FINE)

"""
    _face_categories(mesh, columns, coarse, fault_flux) -> Vector{Int}

The category (0–5, see the file's header) of every face of the `CoarseMesh`
`coarse` built on `mesh` from `columns`; `fault_flux` is `:virtual` or
`:co2lab`. Returns a new vector.
"""
function _face_categories(mesh, columns::VEColumns, coarse, fault_flux::Symbol)
    fault_flux in (:co2lab, :virtual) ||
        throw(ArgumentError("fault_flux: :co2lab or :virtual, got $fault_flux"))
    sealed = Set(columns.sealed)
    top(b) = columns.cells[columns.pos[b]]
    # VE blocks per pillar: co2lab's smooth flux only where a pillar has one
    nve = Dict{Tuple{Int, Int}, Int}()
    for b in eachindex(columns.is_ve)
        ij = columns.column_ij[b]
        columns.is_ve[b] && (nve[ij] = get(nve, ij, 0) + 1)
    end
    single(b) = nve[columns.column_ij[b]] == 1
    map(1:number_of_faces(coarse)) do f
        L, R = coarse.face_neighbors[f]
        F = coarse.coarse_faces_to_fine[f]
        all(in(sealed), F) && return _FACE_SEALED
        vl, vr = columns.is_ve[L], columns.is_ve[R]
        !vl && !vr && return _FACE_FINE
        vl != vr && return _FACE_VE_FINE
        columns.column_ij[L] == columns.column_ij[R] && return _FACE_VE_VERTICAL
        fault_flux == :co2lab && single(L) && single(R) && return _FACE_VE
        # the caprock continues if the blocks' top cells share a face, and the
        # whole-column flux holds if they face each other over their heights
        tops = minmax(top(L), top(R))
        continuous = any(ff -> minmax(mesh.faces.neighbors[ff]...) == tops, F) &&
                     _faces_whole_columns(mesh, columns, F, L, R)
        return continuous ? _FACE_VE : _FACE_VE_DISCONTINUOUS
    end
end

# Whether the fine faces F touch every cell of blocks L and R, i.e. the columns
# face each other over their whole height. If not, the whole-column flux would
# spread a plume's mobility over a height the face does not have (a 20 m
# column beside a 100 m one: its plume flowing as through a fifth of it).
function _faces_whole_columns(mesh, columns, F, L, R)
    touched = Set(c for ff in F for c in mesh.faces.neighbors[ff])
    return all(in(touched), column_cells(columns, L)) &&
           all(in(touched), column_cells(columns, R))
end

"""
    _face_intervals(mesh, columns, coarse, category, geometry) -> (intervals, depth)

The data of each transition face (categories 3–5) for its virtual-cell flux.
Returns a new `4 × nf` matrix and a new vector, `NaN` for other faces:

- `intervals[:, f]`: the depth interval `(t_L, b_L, t_R, b_R)` of each
  side's cells at the face, below that side's block top (`geometry.z_top`).
- `depth[f]`: the common (absolute) depth: the fine cell's centroid for a
  VE–fine face, the area-weighted mean depth of the fine faces otherwise.
"""
function _face_intervals(mesh, columns::VEColumns, coarse, category, geometry)
    nf = number_of_faces(coarse)
    I, depth = fill(NaN, 4, nf), fill(NaN, nf)
    # no transitions (pure VE): no fine geometry needed
    any(in(_TRANSITIONS), category) || return (I, depth)
    geo = tpfv_geometry(mesh)
    zt, zb = _fine_cell_depths(mesh, geo)
    for f in 1:nf
        category[f] in _TRANSITIONS || continue
        F = coarse.coarse_faces_to_fine[f]
        sides = coarse.face_neighbors[f]
        for (s, B) in enumerate(sides)
            # each fine face's cell on this side
            cs = [columns.partition[n[1]] == B ? n[1] : n[2]
                  for n in (mesh.faces.neighbors[ff] for ff in F)]
            I[2s - 1, f] = minimum(zt[cs]) - geometry.z_top[B]
            I[2s, f] = maximum(zb[cs]) - geometry.z_top[B]
        end
        depth[f] = _common_depth(columns, geo, F, sides, category[f])
    end
    return (I, depth)
end

# The depth at which a transition face compares the two sides' potentials
function _common_depth(columns, geo, F, (L, R), category)
    if category == _FACE_VE_FINE
        # the fine cell's centroid (its block has one cell)
        b = columns.is_ve[L] ? R : L
        return geo.cell_centroids[3, columns.cells[columns.pos[b]]]
    end
    a = geo.areas[F]
    return sum(a .* geo.face_centroids[3, F]) / sum(a)
end

"""
    _hybrid_transmissibility(domain, mesh, coarse, category, perm, poro, m)
        -> Union{Vector{Float64}, Nothing}

The coarse faces' transmissibilities by category: 0 for sealed faces, the
coarse mesh's for continuous VE–VE faces (times the area-weighted mean of
their fine faces' multipliers `m`), and for the others the sum of the fine
grid's over the face's fine faces, each times its multiplier (`perm`:
3 × cells of `mesh`).

Returns a new vector, or `nothing` when every face is continuous VE–VE and
every multiplier 1 (the coarse mesh's transmissibilities then apply
unchanged).
"""
function _hybrid_transmissibility(domain, mesh, coarse, category, perm, poro, m)
    plain = all(==(1), m)
    all(==(_FACE_VE), category) && plain && return nothing
    T = JutulDarcy.reservoir_transmissibility(domain)
    if any(c -> c != _FACE_VE && c != _FACE_SEALED, category)
        fine = reservoir_domain(mesh; permeability = perm, porosity = poro)
        Tf = JutulDarcy.reservoir_transmissibility(fine) .* m
    end
    # continuous VE–VE faces: their fine faces' multipliers, weighted by area
    a = plain ? nothing : tpfv_geometry(mesh).areas
    for f in eachindex(category)
        c, F = category[f], coarse.coarse_faces_to_fine[f]
        c == _FACE_SEALED && (T[f] = 0.0)
        c in (_FACE_FINE, _TRANSITIONS...) && (T[f] = sum(Tf[ff] for ff in F))
        c == _FACE_VE && !plain && (T[f] *= sum(a[F] .* m[F]) / sum(a[F]))
    end
    return T
end

# Warn if the blocks fall apart into unconnected pieces (zero transmissibility)
function _warn_disconnected(domain, T)
    N = domain[:neighbors]
    g = SimpleGraph(size(domain[:cell_centroids], 2))
    foreach(f -> T[f] > 0 && add_edge!(g, N[1, f], N[2, f]), axes(N, 2))
    n = length(connected_components(g))
    n > 1 && @warn "The blocks fall apart into $n unconnected pieces; " *
                   "model them separately"
    return nothing
end
