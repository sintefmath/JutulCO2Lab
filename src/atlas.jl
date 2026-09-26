# The Norwegian CO2 Storage Atlas: reading the data sets (roadmap Stage 5,
# step 1) — co2lab-common/atlasgrid (getAtlasGrid.m, helpers/readAAIGrid.m,
# helpers/get*SeaNames.m).
#
# The data sets are MRST's co2atlas, co2atlasnorwegiansea and
# co2atlasbarentssea (NPD data converted to ESRI ASCII grids). The package only
# reads a local copy — it does not download or redistribute them — from
# `atlas_data_directory()`, laid out as MRST's dataset manager expects.
#
# Maps are stored as co2lab stores them: `data[i, j]` with `i` along x (east)
# and `j` along y (north, south first); missing values, infinities and zeros
# are NaN.

# --------------------------------- Data sets ----------------------------------

"The three Atlas data sets, as directory names under the data directory."
const _ATLAS_DATASETS = ("CO2Atlas", "CO2AtlasNorwegianSea", "CO2AtlasBarentsSea")

"""
    atlas_data_directory() -> String

Directory holding the Atlas data sets (`CO2Atlas`, `CO2AtlasNorwegianSea`,
`CO2AtlasBarentsSea`): the environment variable `CO2LAB_ATLAS_DATA` if set,
otherwise `mrst-data` next to the project directory.
"""
atlas_data_directory() =
    get(ENV, "CO2LAB_ATLAS_DATA", normpath(joinpath(@__DIR__, "..", "..", "..",
                                                    "mrst-data")))

"""
    atlas_formation_names(region = :all) -> Vector{String}

co2lab's formation lists (`getNorthSeaNames`, `getNorwegianSeaNames`,
`getBarentsSeaNames`): `region` is `:north_sea`, `:norwegian_sea`,
`:barents_sea` or `:all` (in that order). Returns a new vector.
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

Read an ESRI ASCII grid, as co2lab's `readAAIGrid`: header `ncols`, `nrows`,
`xllcorner`, `yllcorner`, `cellsize`, `NODATA_value`, then the rows north to
south; decimal commas are accepted. NODATA, infinite and zero values become
NaN; rows are flipped to south first and the matrix transposed.

co2lab's reader loses the first data row (its header loop consumes it) and
reads only `nrows − 1` further rows, so its map lies one row north of the
position the header states, with its southernmost row empty.
`co2lab_rows = true` reproduces that, for comparisons with co2lab; in the
Atlas files the first row is always empty, so no values are lost either way.

Returns the NamedTuple `(meta, data)`: `meta = (ncols, nrows, xllcorner,
yllcorner, cellsize, nodata)` as in the header, `data` a new
`Matrix{Float64}` of size `ncols × nrows` (x, then y from the south).
"""
function read_aai_grid(path::AbstractString; co2lab_rows::Bool = false)
    lines = readlines(path)
    header = Dict{String, Float64}()
    k = 0
    while k < length(lines)
        tok = split(replace(lines[k+1], ',' => '.'))
        (length(tok) == 2 && !isnothing(tryparse(Float64, tok[2]))) || break
        header[tok[1]] = parse(Float64, tok[2])
        k += 1
    end

    nx, ny = Int(header["ncols"]), Int(header["nrows"])
    rows = _aai_rows(lines, k, nx, ny, co2lab_rows)

    nodata = header["NODATA_value"]
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

All map files of the three data sets, with co2lab's naming (`readAtlasGrids`
in `getAtlasGrid.m`): the name is the file name without its `top` and
`thickness` parts, joined (`Huginfm_east_thickness` → `Huginfmeast`), and the
variant is `:thickness` if the file name has that part, `:top` otherwise
(this includes rock maps such as `Stofm_poro`, as in co2lab). `.mat`, `.prj`
and README files are skipped.

Returns a new vector of `(name, variant, path)`, data set by data set, files
sorted by name as `dir` lists them.
"""
function atlas_map_files()
    out = NamedTuple{(:name, :variant, :path), Tuple{String, Symbol, String}}[]
    for ds in _ATLAS_DATASETS
        dir = joinpath(atlas_data_directory(), ds)
        isdir(dir) || error("Atlas data set $ds not found in " *
                            "$(atlas_data_directory()) (set CO2LAB_ATLAS_DATA)")
        for f in sort(readdir(dir))
            _is_atlas_map(f) || continue
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
    read_atlas_maps(names = nothing; coarsening = 1, co2lab_rows = false)
        -> Vector{NamedTuple}

The Atlas maps of the formations `names` (all maps for `nothing`), as
co2lab's `readAtlasGrids`: each read with [`read_aai_grid`](@ref); top maps
with any negative value are negated (depth positive down); every
`coarsening`-th value is kept in each direction, from the first, and the cell
size multiplied accordingly. Names compare case-insensitively.

Returns a new vector of `(name, variant, meta, data)`, in the order of
[`atlas_map_files`](@ref); `meta` is the header's with `cellsize` scaled and
`dims = size(data)` added.
"""
function read_atlas_maps(names = nothing; coarsening::Integer = 1,
                         co2lab_rows::Bool = false)
    wanted(n) = isnothing(names) || any(m -> lowercase(m) == lowercase(n), names)
    out = []
    for f in atlas_map_files()
        wanted(f.name) || continue
        g = read_aai_grid(f.path; co2lab_rows)
        # top maps with negative values are elevations: make them depths
        negative = f.variant == :top && any(<(0), filter(!isnan, g.data))
        data = negative ? -g.data : g.data
        data = data[1:coarsening:end, 1:coarsening:end]
        meta = (; g.meta..., cellsize = g.meta.cellsize * coarsening, dims = size(data))
        push!(out, (name = f.name, variant = f.variant, meta = meta, data = data))
    end
    return out
end

# ------------------------------- Average rock ---------------------------------

"""
    atlas_average_rock(name) -> NamedTuple

The Atlas's average rock values for formation `name` (co2lab's `getAvgRock`,
from the Atlas reports; a trailing `fm` is ignored): `(perm, poro, ntg)` with
permeability in m², `NaN` where the Atlas gives none (sealing formations,
unknown names), and `ntg` `nothing` unless the Atlas gives it.
"""
function atlas_average_rock(name::AbstractString)
    n = lowercase(name)
    endswith(n, "fm") && (n = n[1:end-2])
    v = get(_ATLAS_AVERAGE_ROCK, n, (NaN, NaN))
    return (perm = v[1] * 1e-3 * _mrst_darcy(), poro = v[2],
            ntg = length(v) == 3 ? v[3] : nothing)
end

# MRST's darcy (core/utils/units/darcy.m), evaluated in its order so that
# permeabilities agree to the last bit: 1 cP, 1 cm³/s through 1 cm² at 1 atm/cm.
function _mrst_darcy()
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
