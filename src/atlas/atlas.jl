# The Norwegian CO2 Storage Atlas, part 1: reading the data sets (co2lab's
# getAtlasGrid, readAAIGrid, get*SeaNames). The grid is built in
# atlas_grid.jl, the VE domain in atlas_domain.jl.
#
# The data sets are MRST's co2atlas, co2atlasnorwegiansea and
# co2atlasbarentssea (NPD data as ESRI ASCII grids). They are only read from a
# local copy in `atlas_data_directory()`, laid out as MRST's dataset manager
# expects; nothing is downloaded or redistributed.
#
# Map convention, as co2lab: `data[i, j]` with `i` along x (east) and `j`
# along y (south first); missing values, infinities and zeros are NaN.

# --------------------------------- Data sets ----------------------------------

"The three Atlas data sets, as directory names under the data directory."
const _ATLAS_DATASETS = ("CO2Atlas", "CO2AtlasNorwegianSea", "CO2AtlasBarentsSea")

"""
    atlas_data_directory() -> String

Directory holding the Atlas data sets (`CO2Atlas`, `CO2AtlasNorwegianSea`,
`CO2AtlasBarentsSea`), and the Sleipner and Johansen data.

Returns the environment variable `CO2LAB_ATLAS_DATA` if set, otherwise
`mrst-data` next to the project directory.
"""
atlas_data_directory() =
    get(ENV, "CO2LAB_ATLAS_DATA", normpath(joinpath(@__DIR__, "..", "..", "..", "..",
                                                    "mrst-data")))

"""
    atlas_formation_names(region = :all) -> Vector{String}

The Atlas formation names of a region, as co2lab's `getNorthSeaNames`,
`getNorwegianSeaNames` and `getBarentsSeaNames`.

# Arguments

- `region`: `:north_sea`, `:norwegian_sea`, `:barents_sea` or `:all` (the
  three in that order). Any other value throws an `ArgumentError`.

# Returns

A new vector of formation names (e.g. `"Utsirafm"`).
"""
function atlas_formation_names(region::Symbol = :all)
    north = ["Brentgrp", "Brynefm", "Sleipnerfm", "Sognefjordfm", "Fensfjordfm",
             "Krossfjordfm", "Huginfmeast", "Huginfmwest", "Sandnesfm", "Ulafm",
             "Cookfm", "Dunlingp", "Gassumfm", "Johansenfm", "Paleocene",
             "Pliocenesand", "Skadefm", "Statfjordfm", "Utsirafm"]
    norwegian = ["Arefm", "Rorfm", "Notfm", "Ilefm", "Garnfm", "Tiljefm"]
    barents = ["Tubaenfm", "Stofm", "Nordmelafm", "Knurrfm", "Fruholmenfm",
               "Bjarmelandfm"]
    region == :north_sea && return north
    region == :norwegian_sea && return norwegian
    region == :barents_sea && return barents
    region == :all && return vcat(north, norwegian, barents)
    throw(ArgumentError("unknown region :$region"))
end

# ----------------------------- ESRI ASCII grids -------------------------------

"""
    read_aai_grid(path; co2lab_rows = false) -> NamedTuple

Read one ESRI ASCII grid file (an Atlas map), as co2lab's `readAAIGrid`.

The file has a header (`ncols`, `nrows`, `xllcorner`, `yllcorner`,
`cellsize`, `NODATA_value`), then the rows north to south; decimal commas
are accepted.

# Arguments

- `path`: the file.
- `co2lab_rows`: reproduce co2lab's reader, which drops the first data row
  and so places the map one row north, its southernmost row empty (for
  comparisons with co2lab). Default `false`. In the Atlas files the first
  row is always empty, so no values are lost either way.

# Returns

The NamedTuple `(meta, data)`:

- `meta`: the header, `(ncols, nrows, xllcorner, yllcorner, cellsize,
  nodata)`.
- `data`: a new `Matrix{Float64}`, `ncols × nrows` (x, then y from the
  south); NODATA, infinite and zero values are `NaN`.
"""
function read_aai_grid(path::AbstractString; co2lab_rows::Bool = false)
    lines = readlines(path)
    header = Dict{String, Float64}()
    k = 0
    # header lines (keyword, number) until the first data line; some files use
    # decimal commas
    while k < length(lines)
        tok = split(replace(lines[k+1], ',' => '.'))
        (length(tok) == 2 && !isnothing(tryparse(Float64, tok[2]))) || break
        header[tok[1]] = parse(Float64, tok[2])
        k += 1
    end

    nx, ny = Int(header["ncols"]), Int(header["nrows"])
    rows = _aai_rows(lines, k, nx, ny, co2lab_rows)

    nodata = header["NODATA_value"]
    # zero counts as no data too (co2lab's rule)
    data = map(v -> (v == nodata || isinf(v) || v == 0) ? NaN : v, rows)
    meta = (ncols = nx, nrows = ny, xllcorner = header["xllcorner"],
            yllcorner = header["yllcorner"], cellsize = header["cellsize"],
            nodata = nodata)
    # north-first rows → south first, then x along the first index
    return (meta = meta, data = permutedims(data[end:-1:1, :]))
end

# The data rows (ny × nx, north first) after `k` header lines; co2lab's reader
# skips the first data row and leaves the last one empty (zeros).
function _aai_rows(lines, k, nx, ny, co2lab_rows)
    rows = zeros(ny, nx)
    # co2lab skips the first data row, which shifts the map one row north
    first_row = co2lab_rows ? 2 : 1
    for (r, i) in enumerate(first_row:ny)
        vals = parse.(Float64, split(replace(lines[k+i], ',' => '.')))
        rows[r, 1:length(vals)] = vals
    end
    return rows
end

# ------------------------------ Formation maps --------------------------------

"""
    atlas_map_files() -> Vector{NamedTuple}

All map files of the three Atlas data sets, named as co2lab's
`readAtlasGrids` names them. Throws an error if a data set is missing.

# Returns

A new vector of NamedTuples, data set by data set, sorted by file name:

- `name`: the file name without its `top` and `thickness` parts, joined
  (`Huginfm_east_thickness` → `"Huginfmeast"`).
- `variant`: `:thickness` if the file name has that part, otherwise `:top`
  (also for rock maps such as `Stofm_poro`, as in co2lab).
- `path`: the file's path.

`.mat`, `.prj` and README files are skipped.
"""
function atlas_map_files()
    out = NamedTuple{(:name, :variant, :path), Tuple{String, Symbol, String}}[]
    for ds in _ATLAS_DATASETS
        dir = joinpath(atlas_data_directory(), ds)
        isdir(dir) || error("Atlas data set $ds not found in " *
                            "$(atlas_data_directory()) (set CO2LAB_ATLAS_DATA)")
        for f in sort(readdir(dir))
            _is_atlas_map(f) || continue
            # e.g. Huginfm_east_thickness → Huginfmeast, :thickness
            parts = split(f, '_')
            thickness = any(p -> lowercase(p) == "thickness", parts)
            name = join(filter(p -> !(lowercase(p) in ("thickness", "top")), parts))
            push!(out, (name = name, variant = thickness ? :thickness : :top,
                        path = joinpath(dir, f)))
        end
    end
    return out
end

_is_atlas_map(f) = !(lowercase(f) in ("readme.txt", "readme")) &&
                   !(length(f) >= 4 && lowercase(f[end-3:end]) in (".mat", ".prj"))

"""
    read_atlas_maps(names = nothing; coarsening = 1, co2lab_compatible = false)
        -> Vector{NamedTuple}

Read the Atlas maps of some formations, as co2lab's `readAtlasGrids`.
Top maps are returned as depths (positive down): a top map with any
negative value is negated.

# Arguments

- `names`: formation names, compared case-insensitively; `nothing` (the
  default) for all maps.
- `coarsening`: keep every `coarsening`-th value in each direction, from the
  first. Default 1.
- `co2lab_compatible`: place the maps as co2lab does, for comparisons (see
  Notes). Default `false`.

# Returns

A new vector of NamedTuples `(name, variant, meta, data)`, in the order of
[`atlas_map_files`](@ref):

- `name`, `variant`: as [`atlas_map_files`](@ref).
- `meta`: the file header ([`read_aai_grid`](@ref)) with `cellsize` the
  spacing of the values (scaled by `coarsening`), plus `dims = size(data)`
  and `x0`, `y0`, the position (m) of `data[1, 1]`.
- `data`: the values, `Matrix{Float64}` (x, then y from the south).

# Notes

- By default each value sits at the centre of its file cell, so
  `data[1, 1]` is at `(xllcorner, yllcorner)` plus half the file's cell
  size. co2lab instead puts it at the corner `(xllcorner, yllcorner)` and
  reads with its row shift (`co2lab_rows`); `co2lab_compatible = true`
  reproduces both.
"""
function read_atlas_maps(names = nothing; coarsening::Integer = 1,
                         co2lab_compatible::Bool = false)
    wanted(n) = isnothing(names) || any(m -> lowercase(m) == lowercase(n), names)
    out = []
    for f in atlas_map_files()
        wanted(f.name) || continue
        g = read_aai_grid(f.path; co2lab_rows = co2lab_compatible)
        # top maps with negative values are elevations: make them depths
        negative = f.variant == :top && any(<(0), filter(!isnan, g.data))
        data = negative ? -g.data : g.data
        push!(out, (name = f.name, variant = f.variant,
                    _coarsen_map(g.meta, data, coarsening, co2lab_compatible)...))
    end
    return out
end

# A map read by read_aai_grid, coarsened, as (meta, data) with x0, y0 the
# position of its first value (see read_atlas_maps)
function _coarsen_map(meta, data, coarsening, co2lab_compatible)
    data = data[1:coarsening:end, 1:coarsening:end]
    # a value sits at its cell's centre (co2lab: at the corner)
    shift = co2lab_compatible ? 0.0 : meta.cellsize / 2
    meta = (; meta..., cellsize = meta.cellsize * coarsening, dims = size(data),
            x0 = meta.xllcorner + shift, y0 = meta.yllcorner + shift)
    return (meta = meta, data = data)
end

# ------------------------------- Average rock ---------------------------------

"""
    atlas_average_rock(name) -> NamedTuple

The Atlas reports' average rock values for formation `name`, as co2lab's
`getAvgRock`. A trailing `fm` in `name` is ignored, case is not significant.

# Returns

The NamedTuple `(perm, poro, ntg)`:

- `perm`: permeability (m²); `NaN` where the Atlas gives none (sealing
  formations, unknown names).
- `poro`: porosity; `NaN` likewise.
- `ntg`: net-to-gross, or `nothing` unless the Atlas gives it.
"""
function atlas_average_rock(name::AbstractString)
    # the table's names drop the "fm"
    n = lowercase(name)
    endswith(n, "fm") && (n = n[1:end-2])
    v = get(_ATLAS_AVERAGE_ROCK, n, (NaN, NaN))
    # mD to m², with MRST's darcy
    return (perm = v[1] * 1e-3 * _mrst_darcy(), poro = v[2],
            ntg = length(v) == 3 ? v[3] : nothing)
end

# MRST's darcy (m²), evaluated in MRST's order so that permeabilities agree
# to the last bit: 1 cP, 1 cm³/s through 1 cm² at 1 atm/cm.
function _mrst_darcy()
    # computed as MRST's darcy is, for identical rounding
    mu, cm = 1e-2 * (0.1 * 1.0 * 1.0), 1e-2 * 1.0
    p_grad = 101325.0 * 1.0 / cm
    vel = (cm^3 / 1.0) / cm^2
    return vel * mu / p_grad
end

# getAvgRock's table: permeability (mD), porosity, and net-to-gross where given
const _ATLAS_AVERAGE_ROCK = Dict(
    "utsira" => (1000, 0.2112), "skade" => (1000, 0.2112), "sandnes" => (150, 0.0875),
    "sognefjord" => (300, 0.1949), "fensfjord" => (300, 0.1949),
    "krossfjord" => (300, 0.1949), "statfjord" => (200, 0.1071),
    "gassum" => (450, 0.1165), "farsund" => (150, 0.0960), "johansen" => (300, 0.20),
    "cook" => (300, 0.20), "fiskebank" => (1000, 0.2500), "stord" => (15, 0.0600),
    "huginfmeast" => (500, 0.1254), "huginfmwest" => (500, 0.1254),
    "bryne" => (150, 0.1280), "ula" => (300, 0.185), "bjarmeland" => (300, 0.23),
    "sto" => (500, 0.15), "nordmela" => (1, 0.15), "tubaen" => (500, 0.15),
    "tilje" => (140, 0.21, 0.30), "are" => (140, 0.21, 0.30),
    "garn" => (580, 0.27, 0.25), "ile" => (580, 0.27, 0.25))

# ------------------------------ Coastline, wells ------------------------------

"""
    atlas_coastline() -> Vector{Matrix{Float64}}

The coastline of Norway shipped with the North Sea data set (`norway.mat`,
from the Norwegian Mapping Authority), for orientation on Atlas maps.

# Returns

New matrices, one `n × 2` matrix of `(x, y)` (m) per polyline; `NaN` rows
separate pieces. Coordinates are in the Atlas's UTM zone 31.

# Notes

- The file is in UTM zone 33 and is reprojected here. co2lab's
  `showCO2atlas` instead shifts x by 680 km, which misplaces the coast by
  up to tens of kilometres.
"""
function atlas_coastline()
    lines = matread(joinpath(atlas_data_directory(), "CO2Atlas", "norway.mat"))["norway"]
    # the coastline is in UTM zone 33, the Atlas in 31
    rezone(l) = reduce(vcat, [collect(_utm_rezone(l[i, 1], l[i, 2], 33, 31))'
                              for i in axes(l, 1)])
    return [rezone(l) for l in vec(lines)]
end

"""
    atlas_wells() -> NamedTuple

Positions of the wells on the Norwegian continental shelf, shipped with
the North Sea data set (`welldata.mat`, from the NPD fact pages).

Returns the NamedTuple `(x, y)` of new vectors (m, as the file gives them).
"""
function atlas_wells()
    w = matread(joinpath(atlas_data_directory(), "CO2Atlas", "welldata.mat"))["welldata"]
    # stored as (y, x)
    return (x = w[:, 2], y = w[:, 1])
end
