# The Norwegian CO2 Storage Atlas: from top and thickness maps to a corner-point
# grid section (roadmap Stage 5, step 2) — co2lab-common/atlasgrid
# (getAtlasGrid.m, helpers/convertAtlasTo3D.m, helpers/refineInfo.m,
# helpers/updateWithHeterogeneity.m).
#
# The grid's pillars stand at the positions of the top map's values; top and
# thickness (and rock maps) are interpolated bilinearly onto them exactly as
# MATLAB/Octave's `interp2` does, so that co2lab's decks are reproduced to the
# last bit with `co2lab_compatible = true`.

# ----------------------------- Jurassic formations ----------------------------

"Formations whose top is the `Jurassicmid` map (co2lab's `jurassic_mid`)."
const _JURASSIC_MID = ("Brentgrp", "Brynefm", "Sleipnerfm")

"Formations whose top is the `Jurassic` map (co2lab's `jurassic_top`)."
const _JURASSIC_TOP = ("Sognefjordfm", "Fensfjordfm", "Krossfjordfm", "Huginfmeast",
                       "Huginfmwest", "Sandnesfm", "Ulafm")

# ------------------------------ The grid section ------------------------------

"""
    atlas_grid_section(name; coarsening = 1, refining = 1, nz = 1,
                       co2lab_compatible = false) -> NamedTuple

Grid section of Atlas formation `name` (case-insensitive), as co2lab's
`getAtlasGrid(name, 'coarsening', c, 'refining', r, 'nz', nz)`:

- the top map (the formation's own, or the `Jurassic`/`Jurassicmid` map for
  the Jurassic formations) and the thickness map, read with
  [`read_atlas_maps`](@ref) at `coarsening`, then refined `refining` times
  laterally by bilinear interpolation;
- one pillar at each value of the top map; top and thickness interpolated
  bilinearly onto the pillars, the thickness split into `nz` equal layers;
- a cell is active when all four of its pillars have a top and a positive
  thickness;
- rock maps of the formation (`PORO`, `PERMX` = `PERMY` = `PERMZ` in m²,
  `NTG`), where the Atlas has them, interpolated to the cell centres; `NaN`
  where the maps have no value.

`co2lab_compatible = true` reproduces co2lab's deck exactly, for comparisons:
co2lab's placement of the maps ([`read_atlas_maps`](@ref)) and its
activity mask, which also removes the cells east and north of an undefined
pillar and wraps around the lattice's edges (`co2lab-issues.org`).

Returns `(grid, petro)`, both new:
- `grid`: a GeoEnergyIO grid section (`Dict{String, Any}`) for
  `mesh_from_grid_section` — `"cartDims"` `(nx, ny, nz)`, `"COORD"`,
  `"ZCORN"` (both `Vector{Float64}`, 0 at undefined pillars, as GeoEnergyIO's
  `cpgrid_from_horizons`), `"ACTNUM"` (`Array{Bool, 3}`), and the rock arrays
  above (`Array{Float64, 3}`) if present;
- `petro`: the formation's average rock `(perm, poro, ntg)` as
  [`atlas_average_rock`](@ref), with each value replaced by the mean of its
  rock map where there is one (co2lab's `petroinfo`).
"""
function atlas_grid_section(name::AbstractString; coarsening::Integer = 1,
                            refining::Integer = 1, nz::Integer = 1,
                            co2lab_compatible::Bool = false)
    top, thick = _atlas_top_and_thickness(name, coarsening, co2lab_compatible)
    if refining > 1
        top, thick = _refine_map(top, refining), _refine_map(thick, refining)
    end

    Z = _atlas_node_depths(top, thick, nz)
    nx, ny = size(Z, 1) - 1, size(Z, 2) - 1
    active = _atlas_active_cells(Z, co2lab_compatible)
    grid = Dict{String, Any}("cartDims" => (nx, ny, nz),
                             "COORD" => _atlas_coord(top, Z),
                             "ZCORN" => _atlas_zcorn(Z),
                             "ACTNUM" => repeat(active, 1, 1, nz))

    petro = atlas_average_rock(name)
    # co2lab adds rock maps only to formations with their own top map
    if !any(j -> lowercase(j) == lowercase(name), (_JURASSIC_MID..., _JURASSIC_TOP...))
        petro = _add_atlas_rock_maps!(grid, petro, name, top.meta; coarsening,
                                      refining, co2lab_compatible)
    end
    return (grid = grid, petro = petro)
end

# The formation's top and thickness maps (Atlas map NamedTuples, 39 in
# STRUCTURES.md); the Jurassic formations take the Jurassic top maps.
function _atlas_top_and_thickness(name, coarsening, co2lab_compatible)
    is(n) = m -> lowercase(m) == lowercase(n)
    top_name = any(is(name), _JURASSIC_MID) ? "Jurassicmid" :
               any(is(name), _JURASSIC_TOP) ? "Jurassic" : name
    maps = read_atlas_maps([name, top_name]; coarsening, co2lab_compatible)

    find(n, v) = filter(m -> is(n)(m.name) && m.variant == v, maps)
    top, thick = find(top_name, :top), find(name, :thickness)
    (isempty(top) || isempty(thick)) &&
        throw(ArgumentError("no top and thickness maps for formation $name " *
                            "in $(atlas_data_directory())"))
    return (first(top), first(thick))
end

# ------------------------------- Deck arrays ----------------------------------

# Depths of the nz + 1 layer surfaces at each pillar (nx + 1, ny + 1, nz + 1):
# top and thickness / nz interpolated to the top map's positions
# (convertAtlasTo3D); NaN where either is undefined or the thickness is not
# positive.
function _atlas_node_depths(top, thick, nz)
    X, Y = _map_positions(top.meta)
    tx, ty = _map_positions(thick.meta)
    Zt, Zh = permutedims(top.data), permutedims(thick.data)
    Z = zeros(length(X), length(Y), nz + 1)
    for j in eachindex(Y), i in eachindex(X)
        z0 = _interp2(X, Y, Zt, X[i], Y[j])
        dz = _interp2(tx, ty, Zh, X[i], Y[j]) / nz
        dz <= 0 && (dz = NaN)
        # (co2lab then sorts the layers; they are already in order)
        for k in 1:nz+1
            Z[i, j, k] = z0 + dz * (k - 1)
        end
    end
    return Z
end

# COORD: one vertical pillar per node, from the top to the bottom surface
function _atlas_coord(top, Z)
    X, Y = _map_positions(top.meta)
    coord = zeros(6, length(X), length(Y))
    for j in eachindex(Y), i in eachindex(X)
        coord[:, i, j] = [X[i], Y[j], Z[i, j, 1], X[i], Y[j], Z[i, j, end]]
    end
    return _zero_undefined!(vec(coord))
end

# ZCORN: each cell's eight corners from its pillars' surface depths
function _atlas_zcorn(Z)
    # corner index m of 2n → node 1 + m ÷ 2: 1, 2, 2, 3, 3, …, n + 1
    ind(d) = [1 + m ÷ 2 for m in 1:2*(size(Z, d) - 1)]
    return _zero_undefined!(vec(Z[ind(1), ind(2), ind(3)]))
end

"Set the non-finite entries of `v` to 0. Modifies: `v`. Returns `v`."
_zero_undefined!(v) = (v[.!isfinite.(v)] .= 0.0; v)

# Active cells of one layer (nx × ny): all four pillars defined; co2lab's mask
# instead dilates the undefined pillars by one in each direction, wrapping
# around the edges (circshift), and takes cell (i, j) from pillar (i, j).
function _atlas_active_cells(Z, co2lab_compatible)
    undef = isnan.(Z[:, :, 1])
    nx, ny = size(undef, 1) - 1, size(undef, 2) - 1
    if co2lab_compatible
        near = copy(undef)
        for di in -1:1, dj in -1:1
            near .|= circshift(undef, (di, dj))
        end
        return .!near[1:nx, 1:ny]
    end
    return .!(undef[1:nx, 1:ny] .| undef[2:end, 1:ny] .|
              undef[1:nx, 2:end] .| undef[2:end, 2:end])
end

# --------------------------------- Rock maps ----------------------------------

"""
    _add_atlas_rock_maps!(grid, petro, name, top_meta; coarsening, refining,
                          co2lab_compatible) -> NamedTuple

The formation's rock maps (`<name>_poro`, `_perm`, `_ntg`), if any, added to
the grid section as cell arrays, interpolated to the cell centres of the
lattice of `top_meta` (co2lab's `updateWithHeterogeneity`; unlike co2lab,
repeated over all layers). Modifies: `grid` (entries `"PORO"`, `"PERMX"`,
`"PERMY"`, `"PERMZ"`, `"NTG"`). Returns `petro` with each value that has a
map replaced by the map's mean over the lattice's cells (a new NamedTuple).
"""
function _add_atlas_rock_maps!(grid, petro, name, top_meta; coarsening, refining,
                               co2lab_compatible)
    nx, ny, nz = grid["cartDims"]
    targets = Dict("poro" => ("PORO", :poro), "perm" => ("PERMX", :perm),
                   "ntg" => ("NTG", :ntg))
    for f in atlas_map_files()
        parts = split(basename(f.path), '_')
        kind = findfirst(p -> lowercase(p) in ("poro", "perm", "ntg"), parts)
        isnothing(kind) && continue
        lowercase(join(parts[1:end .!= kind])) == lowercase(name) || continue

        g = read_aai_grid(f.path; co2lab_rows = co2lab_compatible)
        m = _coarsen_map(g.meta, g.data, coarsening, co2lab_compatible)
        refining > 1 && (m = _refine_map(m, refining))
        # maps in mD; MRST's convertFrom(data, milli*darcy), before interpolating
        if lowercase(parts[kind]) == "perm"
            m = (; m..., data = m.data .* (1e-3 * _mrst_darcy()))
        end
        layer = _atlas_cell_values(m, top_meta, nx, ny)

        key, field = targets[lowercase(parts[kind])]
        grid[key] = repeat(layer, 1, 1, nz)
        petro = (; petro..., field => _nan_mean(layer))
    end
    if haskey(grid, "PERMX")
        grid["PERMY"], grid["PERMZ"] = copy(grid["PERMX"]), copy(grid["PERMX"])
    end
    return petro
end

# A map interpolated to the cell centres of the lattice of `top_meta` (nx × ny)
function _atlas_cell_values(m, top_meta, nx, ny)
    h, x0, y0 = top_meta.cellsize, top_meta.x0, top_meta.y0
    xc = _linspace(x0 + h / 2, nx * h + x0 - h / 2, nx)
    yc = _linspace(y0 + h / 2, ny * h + y0 - h / 2, ny)
    X, Y = _map_positions(m.meta)
    Z = permutedims(m.data)
    return [_interp2(X, Y, Z, xc[i], yc[j]) for i in 1:nx, j in 1:ny]
end

# MATLAB's mean(v(~isnan(v))) (summation order differs: ~1e-13 relative)
_nan_mean(v) = (w = filter(!isnan, vec(v)); sum(w) / length(w))

# --------------------------- Interpolation helpers ----------------------------

# Positions of a map's values along x and y (co2lab's interpolateData)
function _map_positions(meta)
    nx, ny = meta.dims
    h = meta.cellsize
    return (_linspace(meta.x0, (nx - 1) * h + meta.x0, nx),
            _linspace(meta.y0, (ny - 1) * h + meta.y0, ny))
end

# A map refined r times laterally by bilinear interpolation in index space
# (co2lab's refineInfo); its first value stays in place.
function _refine_map(m, r)
    nx, ny = m.meta.dims
    ix = _linspace(1.0, nx, (nx - 1) * r + 1)
    jy = _linspace(1.0, ny, (ny - 1) * r + 1)
    # interp2(data, Xn, Yn): the rows of `data` (x) are Octave's y
    rows, cols = collect(1.0:nx), collect(1.0:ny)
    data = [_interp2(cols, rows, m.data, jy[q], ix[p])
            for p in eachindex(ix), q in eachindex(jy)]
    meta = (; m.meta..., cellsize = m.meta.cellsize / r, dims = size(data))
    return (; m..., meta = meta, data = data)
end

# Octave's linspace(a, b, n), built symmetrically from both ends
function _linspace(a, b, n)
    n == 1 && return [Float64(b)]
    v = zeros(n)
    v[1], v[n] = a, b
    d = (b - a) / (n - 1)
    for i in 1:(n ÷ 2 - 1)
        v[i+1], v[n-i] = a + i * d, b - i * d
    end
    isodd(n) && (v[n ÷ 2 + 1] = a == -b ? 0.0 : (a + b) / 2)
    return v
end

# Octave's interp2(xs, ys, Z, x, y, 'linear') at one point: Z has rows along
# `ys` and columns along `xs`; NaN outside the data and next to NaN values.
function _interp2(xs, ys, Z, x, y)
    (x < xs[1] || x > xs[end] || y < ys[1] || y > ys[end]) && return NaN
    i = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
    j = clamp(searchsortedlast(ys, y), 1, length(ys) - 1)
    a = Z[j, i]
    b = Z[j, i+1] - a
    c = Z[j+1, i] - a
    d = Z[j+1, i+1] - a - b - c
    tx = (x - xs[i]) / (xs[i+1] - xs[i])
    ty = (y - ys[j]) / (ys[j+1] - ys[j])
    return a + b * tx + c * ty + d * tx * ty
end
