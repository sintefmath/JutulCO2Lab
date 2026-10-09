# The faces of a hybrid VE grid (roadmap Stage 12 items 1–2; design:
# analysis/stage12-transition-flux-design.md §3.3).
#
# Every coarse face gets a category:
#   0 sealed — all its fine faces sealed: no flow;
#   1 fine–fine — the ordinary 3D flux;
#   2 VE–VE, continuous — the two blocks' top cells share a face (the caprock
#     continues): the ordinary whole-column VE flux;
#   3 VE–VE, discontinuous — lateral, top cells not connected (a fault offset
#     larger than the top cell, a pinch-out, stacked blocks): the virtual cell;
#   4 VE–VE, vertical — stacked blocks through a leaky layer: the virtual cell
#     over the blocks' adjacent cells;
#   5 VE–fine, lateral or vertical — the virtual cell over the fine cell's
#     extent.
# Categories 3–5 (transitions) carry, per side, the depth interval of that
# side's cells at the face (below its block's top) and a common depth (the
# fine cell's centroid for 5, the faces' mean depth for 3 and 4), and the
# transmissibility summed from the fine grid's. Category 2 keeps the coarse
# mesh's transmissibility (as the pure VE gates), 1 takes the fine grid's.
#
# `fault_flux = :virtual` (default; researcher, 2026-10-05) follows the
# topology; `:co2lab` makes lateral VE–VE faces between pillars with one VE
# block each continuous (co2lab's smooth flux across faults, for the gates
# against MRST; stacked blocks follow the topology in any case).

const _FACE_SEALED = 0
const _FACE_FINE = 1
const _FACE_VE = 2
const _FACE_VE_DISCONTINUOUS = 3
const _FACE_VE_VERTICAL = 4
const _FACE_VE_FINE = 5
const _TRANSITIONS = (_FACE_VE_DISCONTINUOUS, _FACE_VE_VERTICAL, _FACE_VE_FINE)

"""
    _face_categories(mesh, columns, coarse, fault_flux) -> Vector{Int}

The category (0–5, see the file's header) of every face of the `CoarseMesh`
`coarse` built on `mesh` from `columns`. Returns a new vector; the arguments
are only read.
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
        # the caprock continues if the blocks' top cells share a face
        tops = minmax(top(L), top(R))
        continuous = any(ff -> minmax(mesh.faces.neighbors[ff]...) == tops, F)
        return continuous ? _FACE_VE : _FACE_VE_DISCONTINUOUS
    end
end

"""
    _face_intervals(mesh, columns, coarse, category, geometry) -> (intervals, depth)

For each transition face (categories 3–5): `intervals[:, f]` the depth
interval `(t_L, b_L, t_R, b_R)` of each side's cells at the face, below that
side's block top (`geometry.z_top`), and `depth[f]` the common depth
(absolute): the fine cell's centroid for a VE–fine face, the area-weighted
mean depth of the fine faces otherwise. `NaN` for other faces. Returns a new
`4 × nf` matrix and a new vector; the arguments are only read.
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
    _hybrid_transmissibility(domain, mesh, coarse, category, perm, poro, ntg)
        -> Union{Vector{Float64}, Nothing}

The coarse faces' transmissibilities by category: 0 for sealed faces, the
coarse mesh's for continuous VE–VE faces, the sum of the fine grid's over
the face's fine faces for the others (`perm` 3 × cells of `mesh`). `nothing`
when every face is continuous VE–VE (the coarse mesh's apply unchanged).
Returns a new vector; the arguments are only read.
"""
function _hybrid_transmissibility(domain, mesh, coarse, category, perm, poro, ntg)
    all(==(_FACE_VE), category) && return nothing
    T = JutulDarcy.reservoir_transmissibility(domain)
    if any(c -> c != _FACE_VE && c != _FACE_SEALED, category)
        kw = isnothing(ntg) ? (;) : (net_to_gross = ntg,)
        fine = reservoir_domain(mesh; permeability = perm, porosity = poro, kw...)
        Tf = JutulDarcy.reservoir_transmissibility(fine)
    end
    for f in eachindex(category)
        c = category[f]
        c == _FACE_SEALED && (T[f] = 0.0)
        c in (_FACE_FINE, _TRANSITIONS...) &&
            (T[f] = sum(Tf[ff] for ff in coarse.coarse_faces_to_fine[f]))
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
