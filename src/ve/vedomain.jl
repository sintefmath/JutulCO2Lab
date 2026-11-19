# The VE reservoir domain: a CoarseMesh with one block per column, carrying the
# column-averaged rock and the data the VE model needs; also the initial
# pressure and the domain's outer lateral boundary sides. Everything is derived
# from the 3D mesh and its VEColumns partition.

# ------------------------- Rock and initial pressure --------------------------

"""
    net_rock(mesh; poro, perm_x, perm_y = nothing, perm_z = nothing, ntg,
             ntg_min = 0.0) -> NamedTuple

Fold net-to-gross into the rock of a 3D mesh and remove the cells without net
rock.

# Arguments

- `mesh`: the 3D mesh, with logical indices (as for [`ve_columns`](@ref)).
- `poro`: porosity per cell, of the whole (gross) cell.
- `perm_x`, `perm_y`: lateral permeability (m²) per cell; `perm_y` optional.
- `perm_z`: vertical permeability (m²) per cell; optional.
- `ntg`: net-to-gross per cell (≥ 0), or `nothing` (the rock as given, no
  cell removed).
- `ntg_min`: cells with `ntg ≤ ntg_min` are removed. Default 0.

# Returns

The NamedTuple `(mesh, cells, poro, perm_x, perm_y, perm_z)`:

- `mesh`: the mesh without the removed cells ([`remove_cells`](@ref)); `mesh`
  itself, not a copy, when no cell is removed.
- `cells`: the kept cells, as indices into the given `mesh` (for other
  per-cell data).
- `poro`, `perm_x`, `perm_y`: new vectors per kept cell, `poro · ntg`,
  `perm_x · ntg`, `perm_y · ntg`; `perm_y` is `nothing` if not given.
- `perm_z`: a new vector per kept cell, the given `perm_z` **not** scaled;
  `nothing` if not given.

Throws an `ArgumentError` if `ntg` has a wrong length, a negative or a `NaN`
value.

# Notes

- Net-to-gross is the fraction of a cell's thickness that is reservoir rock;
  the rest holds no pore space and lets nothing through laterally. Pore
  volume and lateral flow scale with it, as with ECLIPSE's `NTG`.
- Vertical permeability: a given `perm_z` is taken as effective already, the
  non-reservoir part accounted for (ECLIPSE's `NTG` leaves vertical
  transmissibility alone). Without it, [`ve_domain`](@ref) takes vertical
  permeability from `perm_x`, so one permeability per cell is scaled in all
  directions. JutulDarcy's `net_to_gross` scales lateral transmissibility
  only, also then.
- A removed cell is a hole, with no flow through it in any direction:
  [`ve_columns`](@ref) cuts the pillar there (by default a column above it
  and one below).
- co2lab averages porosity and net-to-gross over a column separately and
  multiplies the averages; that misses the pore volume where net-to-gross
  varies within a column. The average of the product, as here, is exact.
"""
function net_rock(mesh; poro, perm_x, perm_y = nothing, perm_z = nothing, ntg,
                  ntg_min::Real = 0.0)
    nc = number_of_cells(mesh)
    # no net-to-gross: every cell kept, the rock as given
    isnothing(ntg) && return (mesh = mesh, cells = collect(1:nc), poro = copy(poro),
                              perm_x = copy(perm_x), perm_y = _copy_or(perm_y),
                              perm_z = _copy_or(perm_z))
    length(ntg) == nc || throw(ArgumentError("ntg: one value per cell expected"))
    all(>=(0), ntg) || throw(ArgumentError("ntg must be ≥ 0 and not NaN"))

    # cells with net rock kept, the others removed from the mesh
    cells = findall(>(ntg_min), ntg)
    trimmed = remove_cells(mesh, setdiff(1:nc, cells))

    # pore volume and lateral flow scale with net-to-gross; a given vertical
    # permeability is effective already
    net(x) = isnothing(x) ? nothing : x[cells] .* ntg[cells]
    return (mesh = trimmed, cells = cells, poro = net(poro), perm_x = net(perm_x),
            perm_y = net(perm_y), perm_z = isnothing(perm_z) ? nothing : perm_z[cells])
end

# A copy of x, or nothing
_copy_or(x) = isnothing(x) ? nothing : copy(x)

"""
    ve_average_rock(columns, dz; perm_x, perm_y = nothing, poro,
                    anisotropy_rtol = 1e-6) -> NamedTuple

Thickness-weighted column averages of the rock, `sum(x .* dz) / sum(dz)` per
column, exactly as co2lab's `averageRock`. Used by [`ve_domain`](@ref).

# Arguments

- `columns`: the column partition ([`ve_columns`](@ref)).
- `dz`: thickness (m) per fine cell, in `columns.cells` order
  ([`ve_column_geometry`](@ref)).
- `perm_x`, `perm_y`: lateral permeability (m²) per fine cell of the 3D mesh;
  `perm_y` is only checked against `perm_x`.
- `poro`: porosity per fine cell.
- `anisotropy_rtol`: relative tolerance for `perm_y` differing from `perm_x`.

# Returns

The NamedTuple `(perm, poro)` of new per-column vectors.

# Notes

- VE permeability is a **scalar from lateral permeability only**: the
  thickness-weighted x-permeability, as in co2lab. Vertical permeability is
  never used.
- If `perm_y` differs from `perm_x` beyond `anisotropy_rtol` in any cell, a
  warning reports how many cells and the largest ratio; x is used anyway.
- Net-to-gross: average the net rock ([`net_rock`](@ref)).
"""
function ve_average_rock(columns::VEColumns, dz::AbstractVector;
                         perm_x, perm_y = nothing, poro, anisotropy_rtol = 1e-6)
    isnothing(perm_y) || _warn_lateral_anisotropy(perm_x, perm_y, anisotropy_rtol)

    # thickness-weighted average over each column's cells
    avg(x) = [sum(x[columns.cells[r]] .* dz[r]) / sum(dz[r])
              for r in (columns.pos[c]:(columns.pos[c+1] - 1)
                        for c in 1:number_of_columns(columns))]
    return (perm = avg(perm_x), poro = avg(poro))
end

"""
    ve_hydrostatic_pressure(geometry; rho_w, g = Jutul.gravity_constant)
        -> Vector{Float64}

Initial VE pressure per column: the hydrostatic brine pressure at the
**caprock** for a constant brine density, as co2lab's examples initialise
(no surface pressure offset, no JutulDarcy equilibration).

# Arguments

- `geometry`: column geometry ([`ve_column_geometry`](@ref)); uses `z_top`,
  the caprock depth (m, positive down).
- `rho_w`: brine density (kg/m³).
- `g`: gravity (m/s²).

# Returns

A new `Vector{Float64}`, `rho_w · g · z_top` (Pa), one value per column. It
is the brine (reference phase) pressure; with no CO2 present the CO2 pressure
equals it.

# Notes

- Do not compute it from the domain's cell centroids: they are mid-column on
  the `CoarseMesh` and would give a wrong pressure without any error.
- For a pressure-dependent brine density use
  [`ve_hydrostatic_profile`](@ref).
"""
ve_hydrostatic_pressure(geometry; rho_w, g = Jutul.gravity_constant) =
    rho_w * g .* geometry.z_top

"""
    ve_hydrostatic_profile(fluid, ve; g = Jutul.gravity_constant, dz = 1.0)
        -> Function

Hydrostatic brine pressure as a function of depth, consistent with the
fluid's own (possibly pressure- and temperature-dependent) brine density.
Brine at these pressures — at the caprocks, and on open sides at their
depths — is at rest, which a constant `ρ_w` other than the fluid's would not
give. With a constant density it equals [`ve_hydrostatic_pressure`](@ref).

# Arguments

- `fluid`: the fluid ([`co2lab_ve_fluid`](@ref)); its brine density is used.
- `ve`: the formation ([`ve_domain`](@ref)); its caprock depths and column
  temperatures are used.
- `g`: gravity (m/s²).
- `dz`: integration step (m).

# Returns

A new function `z -> p` (depth in m, positive down; pressure in Pa): linear
between the steps, extrapolated with the last step's gradient below them.

# Notes

- Integrates `dp/dz = g ρ_w(p, T(z))` from `p = 0` at `z = 0` (co2lab's
  convention, no surface pressure) with Heun's method, down to below the
  deepest column.
- `T(z)` interpolates the columns' temperatures linearly in caprock depth,
  constant beyond the shallowest and deepest; JutulDarcy's default 30 °C if
  the domain has none.
"""
function ve_hydrostatic_profile(fluid, ve; g = Jutul.gravity_constant, dz = 1.0)
    z_top, d = ve.geometry.z_top, ve.domain
    T = haskey(d, :temperature, Cells()) ? d[:temperature] :
        fill(Jutul.default_value(nothing, JutulDarcy.Temperature()), length(z_top))
    Tz = _depth_interpolant(z_top, T)
    rho(p, z) = _phase_value(fluid.rho_ref[1], fluid.rho_pvt[1], p, Tz(z))

    # Heun's method from the surface to below the deepest caprock
    n = max(ceil(Int, (maximum(z_top) + maximum(ve.geometry.H)) / dz), 1)
    zs, ps = collect(range(0.0, n * dz; length = n + 1)), zeros(n + 1)
    for k in 1:n
        r = rho(ps[k], zs[k])
        ps[k+1] = ps[k] + g * dz * (r + rho(ps[k] + g * dz * r, zs[k+1])) / 2
    end
    return _depth_interpolant(zs, ps; extrapolate = true)
end

# A new function z -> value, linear between the points (z, v) (values at equal
# z averaged); beyond them constant, or continued linearly if `extrapolate`
function _depth_interpolant(z, v; extrapolate = false)
    zu, vu, nu = Float64[], Float64[], Int[]
    for i in sortperm(z)
        if isempty(zu) || z[i] != zu[end]
            push!(zu, z[i]); push!(vu, 0.0); push!(nu, 0)
        end
        vu[end] += v[i]; nu[end] += 1
    end
    vu ./= nu
    length(zu) == 1 && return _ -> vu[1]
    function f(x)
        k = clamp(searchsortedlast(zu, x), 1, length(zu) - 1)
        x = extrapolate ? x : clamp(x, zu[1], zu[end])
        w = (x - zu[k]) / (zu[k+1] - zu[k])
        return (1 - w) * vu[k] + w * vu[k+1]
    end
    return f
end

# Warn if kx and ky differ beyond rtol in any cell; `what` (VE permeability,
# or the integrated relperm) then uses kx (see ve_average_rock)
function _warn_lateral_anisotropy(kx, ky, rtol; what = "VE permeability")
    rel = abs.(kx .- ky) ./ max.(abs.(kx), abs.(ky), floatmin())
    n = count(>(rtol), rel)
    n == 0 && return nothing

    ratio = maximum(max.(kx ./ ky, ky ./ kx))
    @warn "Lateral permeability is anisotropic in $n of $(length(kx)) cells " *
          "(largest kx/ky or ky/kx: $(round(ratio, sigdigits = 3))); $what " *
          "uses x-permeability only"
    return nothing
end

# ------------------------------- The VE domain --------------------------------

"""
    ve_domain(mesh, columns; perm_x, perm_y = nothing, perm_z = nothing, poro,
              rugosity = 0.0, temperature = nothing, fault_flux = :virtual)
        -> NamedTuple

Build the VE reservoir domain: a `CoarseMesh` over `mesh` with one block per
column of `columns`, carrying the column-averaged rock and the data the VE
model needs. Its `domain` is what [`setup_co2lab_model`](@ref) takes.

`columns` may hold stacked columns (`ve_columns(...; below_holes = :stack)`)
and fine 3D cells (hybrid grids); both are handled in full.

# Arguments

- `mesh`, `columns`: the 3D mesh and its column partition
  ([`ve_columns`](@ref)).
- `perm_x`, `poro`: lateral permeability (m²) and porosity per fine cell,
  averaged over each column by thickness.
- `perm_y`: only checked; VE uses `perm_x` (a warning if they differ).
- `perm_z`: vertical permeability per fine cell, for vertical faces between
  blocks (stacked blocks, fine cells); default `perm_x`.
- `rugosity`: the caprock's rugosity `dh` (m), a number or one per column;
  the relperm loses an accretion layer `min(dh, H)` under the caprock
  (co2lab's `makeVEFluid` option).
- `temperature`: the formation's temperature (K) at each column's caprock,
  for tabulated fluid properties: a number, one per column, or a function of
  depth (e.g. [`geothermal_gradient`](@ref)). Default: the model's, 30 °C.
- `fault_flux`: the flux across a face where two columns' caprocks do not
  meet (a fault): `:virtual` compares the two sides at a common depth;
  `:co2lab` joins them as one continuous caprock, as co2lab does, between
  pillars with one block each (for comparisons with co2lab).

# Returns

The NamedTuple `(domain, coarse, geometry, rock)`, all new:

- `domain`: the `DataDomain` on the `CoarseMesh`, with
  - per block: the averaged `:permeability` and `:porosity`;
    `:column_height`, `:rugosity`, `:temperature` (if given);
    `:is_ve`; `:reference_depth`, the depth of the block's pressure (a VE
    block's caprock, a fine cell's centroid); `:pillar_columns`, the number
    of blocks in its pillar;
  - per face: `:gdz`, the gravity term from the reference depths;
    `:face_category`, `:face_intervals`, `:face_depth` (below);
    `:transmissibilities` if some face is not a continuous VE–VE face or
    some multiplier is not 1.
- `coarse`: the `CoarseMesh` that `domain` wraps.
- `geometry`: caprock depth, height and extent per column
  ([`ve_column_geometry`](@ref)).
- `rock`: the averaged rock ([`ve_average_rock`](@ref)).

# Faces

Each face has a category: 0 sealed, 1 fine–fine, 2 continuous VE–VE, and
the transitions 3 discontinuous VE–VE (the columns overlap over part of
their heights: stacked columns, faults), 4 vertical VE–VE (through a leaky
layer), 5 VE–fine. A transition face carries each side's depth interval
(`:face_intervals`) and a common depth (`:face_depth`) for its flux; other
faces have `NaN` there.

Transmissibilities: 0 across sealed faces, the coarse mesh's for continuous
VE–VE faces, otherwise the sum over the fine faces. The columns'
transmissibility multipliers ([`ve_columns`](@ref)) apply per fine face (on
continuous VE–VE faces, their area-weighted mean). A warning is given if
the blocks fall apart into unconnected pieces.

# Notes

- Permeability is not multiplied by the column height: the coarse mesh's
  face areas already contain it, so Jutul's transmissibility equals MRST's
  `perm·H` one.
- Block centroids stay mid-column, as the transmissibilities are computed
  from them. The caprock enters through `:gdz` and `:reference_depth`;
  anything needing a block's pressure depth must use `:reference_depth`,
  not the centroid.
- Net-to-gross: pass the net rock of [`net_rock`](@ref) (with the mesh it
  returns), so the domain carries no `:net_to_gross`.
"""
function ve_domain(mesh, columns::VEColumns; perm_x, perm_y = nothing,
                   perm_z = nothing, poro, rugosity = 0.0,
                   temperature = nothing, fault_flux::Symbol = :virtual)
    geometry = ve_column_geometry(mesh, columns)
    rock = ve_average_rock(columns, geometry.dz; perm_x, perm_y, poro)

    # the VE grid: one coarse cell per column
    coarse = CoarseMesh(mesh, columns.partition)
    domain = reservoir_domain(coarse; permeability = rock.perm, porosity = rock.poro)

    # gravity between VE blocks from caprock depths, not mid-column centroids;
    # fine cells from their centroids
    zref = [columns.is_ve[c] ? geometry.z_top[c] : geometry.z_cell[columns.pos[c]]
            for c in 1:number_of_columns(columns)]
    domain[:gdz, Faces()] = Jutul.compute_face_gdz(domain[:neighbors], zref)
    domain[:reference_depth, Cells()] = zref
    domain[:column_height, Cells()] = geometry.H
    ncol = number_of_columns(columns)
    domain[:rugosity, Cells()] = rugosity isa Number ? fill(Float64(rugosity), ncol) :
                                 collect(Float64, rugosity)
    domain[:is_ve, Cells()] = copy(columns.is_ve)
    domain[:pillar_columns, Cells()] = _pillar_columns(columns)
    isnothing(temperature) ||
        (domain[:temperature, Cells()] = _column_temperature(temperature, geometry))

    # the faces: categories, virtual-cell intervals, transmissibilities
    perm = [perm_x'; something(perm_y, perm_x)'; something(perm_z, perm_x)']
    _set_face_data!(domain, mesh, columns, coarse, geometry, perm, poro, fault_flux)
    return (domain = domain, coarse = coarse, geometry = geometry, rock = rock)
end

# Temperature per column: a scalar, one value per column, or a function of
# the caprock depth. Returns a new vector.
_column_temperature(T::Number, geometry) = fill(Float64(T), length(geometry.H))
_column_temperature(T::AbstractVector, geometry) = collect(Float64, T)
_column_temperature(T::Function, geometry) = Float64.(T.(geometry.z_top))

"""
    geothermal_gradient(; T_ref, z_ref = 0.0, gradient) -> Function

A linear geothermal temperature profile, for [`ve_domain`](@ref)'s
`temperature`.

# Arguments

- `T_ref`: temperature (K) at depth `z_ref`.
- `z_ref`: reference depth (m, positive down).
- `gradient`: temperature gradient (K/m), e.g. 0.03 for 30 °C/km.

# Returns

A new function `z -> T_ref + gradient · (z − z_ref)` of depth `z` (m,
positive down), giving temperature in K.
"""
geothermal_gradient(; T_ref, z_ref = 0.0, gradient) =
    z -> T_ref + gradient * (z - z_ref)

# Number of columns in each column's pillar.
function _pillar_columns(columns::VEColumns)
    n = Dict{Tuple{Int, Int}, Int}()
    foreach(p -> n[p] = get(n, p, 0) + 1, columns.column_ij)
    return [n[p] for p in columns.column_ij]
end

"""
    _set_face_data!(domain, mesh, columns, coarse, geometry, perm, poro, fault_flux)

Store each coarse face's category, virtual-cell intervals, common depth and
transmissibility (see `vehybrid.jl`), and warn if the blocks fall apart.
Returns `domain`.

Modifies: `domain` (sets `:face_category`, `:face_intervals`, `:face_depth`
on `Faces()`, and `:transmissibilities` unless every face is continuous
VE–VE).
"""
function _set_face_data!(domain, mesh, columns, coarse, geometry, perm, poro,
                         fault_flux)
    category = _face_categories(mesh, columns, coarse, fault_flux)
    intervals, depth = _face_intervals(mesh, columns, coarse, category, geometry)
    domain[:face_category, Faces()] = category
    domain[:face_intervals, Faces()] = intervals
    domain[:face_depth, Faces()] = depth
    T = _hybrid_transmissibility(domain, mesh, coarse, category, perm, poro,
                                 columns.multiplier)
    isnothing(T) || (domain[:transmissibilities, Faces()] = T)
    # stacked or sealed blocks may fall apart
    stacked = any(>(1), domain[:pillar_columns])
    if stacked || !isnothing(T)
        Tw = something(T, JutulDarcy.reservoir_transmissibility(domain))
        _warn_disconnected(domain, Tw)
    end
    return domain
end

# ------------------------------- Boundary sides -------------------------------

"""
    VEBoundarySides

The outer lateral boundary sides of the VE columns, one entry per (column,
side) — the counterpart of the boundary edges of co2lab's top-surface grid.
Made by [`ve_boundary_sides`](@ref); boundary conditions are set on them
([`ve_boundary_conditions`](@ref)).

# Fields

- `column`: the column (coarse cell) the side belongs to.
- `direction`: 1..4 = `:left`, `:right`, `:lower`, `:upper` (Jutul's
  direction tags: logical i−, i+, j−, j+).
- `edge_depth`: caprock depth (m) along the side, the mean depth of the two
  top-face nodes bounding it (co2lab's `Gt.faces.z`, without co2lab's
  averaging of nodes across pillar offsets). For a fine cell of a hybrid
  grid: the depth of the side's centroid.
- `edge_midpoint`: (x, y) midpoint of that top edge (fine cell: of the side's
  centroid).
- `area`: the side's area (m²), the sum of its fine boundary faces' areas.
- `centroid`: the side's centroid (m), area-weighted over its fine faces.
- `normal`: the side's outward unit normal, area-weighted over its fine
  faces.
- `half_trans`: two-point half transmissibility (m³) from the column's
  centroid to the side, `k·A·(c·n)/|c|²`, the same formula family as the
  `CoarseMesh`'s interior transmissibilities.
- `fine_faces`: indices of the side's boundary faces in the 3D mesh.
- `across_fault`: a column exists at the logical neighbour position (at any
  depth), but no fine face connects to it: a fault with complete offset. Such
  a side is an outer side, as in co2lab's top-surface grid, but physically
  sealed: give it a boundary condition only deliberately (co2lab's Johansen
  example keeps its main fault sealed).
"""
struct VEBoundarySides
    column::Vector{Int}                        # the column (coarse cell)
    direction::Vector{Int}                     # 1..4: left, right, lower, upper
    edge_depth::Vector{Float64}                # caprock (fine: centroid) depth of side
    edge_midpoint::Vector{NTuple{2, Float64}}  # (x, y) midpoint of the top edge
    area::Vector{Float64}                      # side area (m²)
    centroid::Vector{NTuple{3, Float64}}       # side centroid
    normal::Vector{NTuple{3, Float64}}         # outward unit normal
    half_trans::Vector{Float64}                # half transmissibility to the side
    fine_faces::Vector{Vector{Int}}            # boundary-face indices in the 3D mesh
    across_fault::Vector{Bool}                 # sealed full-offset fault side
end

Base.length(s::VEBoundarySides) = length(s.column)

"""
    ve_boundary_sides(mesh, columns, domain) -> (sides, unassigned)

The outer boundary sides of the columns: the column sides
([`ve_column_sides`](@ref)) with no neighbouring column, where open
boundary conditions can be placed. Which sides are outer is decided by
topology alone (face directions and cell connectivity); geometry only
supplies the quantities.

# Arguments

- `mesh`, `columns`: the 3D mesh and its column partition
  ([`ve_columns`](@ref)).
- `domain`: the VE `DataDomain` ([`ve_domain`](@ref)); its permeability and
  block centroids give the half transmissibilities.

# Returns

The tuple `(sides, unassigned)`:

- `sides`: a new [`VEBoundarySides`](@ref), one entry per outer side, ordered
  by column and, within a column, by direction (left, right, lower, upper).
- `unassigned`: a new `Vector{Int}` of boundary-face indices of `mesh` on
  sides that *do* have a neighbour, in the same order. They stay no-flow.

# Notes

- A column that has a neighbour on some side but stands higher (an internal
  caprock step) exposes boundary faces of the 3D grid there, yet the side is
  not outer: co2lab's top-surface grid has an interior edge there. Such faces
  go to `unassigned`, no-flow as in co2lab. Likewise a stacked column's side
  at the depth of a shale missing from the grid.
- `across_fault` is set where the logical neighbour pillar holds any column,
  at any depth: across a fault with complete offset the columns usually do
  not overlap in depth.
- A side's top edge is where the column's top face meets the top cell's
  boundary faces on that side; it must have exactly two nodes (an error is
  thrown otherwise).
"""
function ve_boundary_sides(mesh::UnstructuredMesh, columns::VEColumns, domain)
    cs = ve_column_sides(mesh, columns)
    g = _boundary_geometry(mesh)
    bdir = _boundary_directions(mesh)
    ijk = [cell_ijk(mesh, c) for c in 1:number_of_cells(mesh)]
    perm, cc = domain[:permeability], domain[:cell_centroids]
    out = (column = Int[], dir = Int[], depth = Float64[], mid = NTuple{2, Float64}[],
           area = Float64[], cen = NTuple{3, Float64}[], nrm = NTuple{3, Float64}[],
           T = Float64[], faces = Vector{Int}[], fault = Bool[])

    # pillars that hold a column (for the fault check below)
    occupied = Set(columns.column_ij)

    for s in 1:length(cs)
        # outer sides only: no column connected beyond
        cs.neighbour[s] == 0 || continue
        col, fs = cs.column[s], cs.boundary_faces[s]

        d, top = cs.direction[s], first(column_cells(columns, col))
        centroid = (cc[1, col], cc[2, col], cc[3, col])
        A, C, N, t = _side_aggregate(g, fs, perm[col], centroid)
        # a fine cell's side: its faces' centroid, the depth its flux uses
        depth, mid = columns.is_ve[col] ?
                     _side_top_edge(g, bdir, _top_face_nodes(mesh, ijk, bdir, top), top,
                                    d, col) : (C[3], (C[1], C[2]))

        push!(out.column, col); push!(out.dir, d)
        push!(out.depth, depth); push!(out.mid, mid)
        push!(out.area, A); push!(out.cen, C); push!(out.nrm, N)
        push!(out.T, t); push!(out.faces, fs)

        # a column exists across this side, but nothing connects to it: sealed fault
        across = _logical_neighbour_pillar(columns.column_ij[col], d)
        push!(out.fault, across in occupied)
    end

    sides = VEBoundarySides(out.column, out.dir, out.depth, out.mid, out.area, out.cen,
                            out.nrm, out.T, out.faces, out.fault)
    return (sides, copy(cs.exposed_faces))
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
    _side_top_edge(g, bdir, top_nodes, top, d, col) -> (depth, midpoint)

Caprock edge of column `col`'s side in direction `d`: the two nodes shared by
the top cell `top`'s top face (`top_nodes`: a boundary face, or for a column
under another block the interior face to the cell above) and its boundary
faces on that side. Returns their mean depth and (x, y) midpoint; throws if
there are not exactly two such nodes.
"""
function _side_top_edge(g, bdir, top_nodes, top, d, col)
    side_nodes = Set{Int}()
    for f in g.c2f[top]
        bdir[f] == d && union!(side_nodes, g.nodes[f])
    end

    shared = [n for n in top_nodes if n in side_nodes]
    length(shared) == 2 ||
        error("column $col, $(_DIRECTIONS[d]) side: top edge has " *
              "$(length(shared)) nodes, expected 2")

    Xa, Xb = g.X[shared[1]], g.X[shared[2]]
    return ((Xa[3] + Xb[3]) / 2, ((Xa[1] + Xb[1]) / 2, (Xa[2] + Xb[2]) / 2))
end

# One side's fine faces fs aggregated: (area, centroid, outward unit normal,
# two-point half transmissibility k·A·(c·n)/|c|² from the column centroid)
function _side_aggregate(g, fs, perm, col_centroid)
    A = sum(g.ba[f] for f in fs)
    # area-weighted centroid and normal of the fine faces
    C = ntuple(d -> sum(g.ba[f] * g.bc[d, f] for f in fs) / A, 3)
    Nv = ntuple(d -> sum(g.ba[f] * g.bn[d, f] for f in fs), 3)
    N = Nv ./ sqrt(sum(abs2, Nv))

    cvec = C .- col_centroid
    t = perm * A * sum(cvec .* N) / sum(abs2, cvec)
    return (A, C, N, t)
end
