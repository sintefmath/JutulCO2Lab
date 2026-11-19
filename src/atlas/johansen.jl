# The Johansen formation's NPD models, both faulted corner-point grids:
#
# - the sector model (NPD5) and the full-field model as meshes of the Johansen
#   layers, as co2lab's trapsJohansen example reads them (for trap analysis);
# - the sector model as a VE domain with its rock, as co2lab's
#   exampleVEJohansen builds it (makeJohansenVEgrid).

# ---------------------------------- The data ----------------------------------

"Deck files and Johansen layers of the two NPD models."
const _JOHANSEN_MODELS = Dict(
    :sector => (file = "NPD5.grdecl", layers = 6:11,
                url = "https://www.sintef.no/project/MatMoRA/Johansen/NPD5.zip"),
    :full_field => (file = "FULLFIELD_IMAXJMAX.GRDECL", layers = 10:14,
                    url = "https://www.sintef.no/project/MatMoRA/Johansen/" *
                          "FULLFIELD_Eclipse.zip"))

"""
    johansen_data_directory() -> String

Directory holding the Johansen NPD decks (`NPD5.grdecl`,
`FULLFIELD_IMAXJMAX.GRDECL` and its fault include file) and the sector
model's rock (`NPD5_Permeability.txt`, `NPD5_Porosity.txt`).

Returns `Johansen` in [`atlas_data_directory`](@ref), as MRST's data
directory has it.
"""
johansen_data_directory() = joinpath(atlas_data_directory(), "Johansen")

# ------------------------------- The grids ------------------------------------

"""
    johansen_mesh(model = :sector; co2lab_compatible = false) -> NamedTuple

The Johansen formation's layers in one of the NPD models, as a mesh with
columns (no rock), cut from the deck as co2lab's `trapsJohansen` cuts them
(`cutGrdecl`). For trap analysis and maps.

# Arguments

- `model`: `:sector` (NPD5, deck layers 6–11) or `:full_field` (deck layers
  10–14); otherwise an `ArgumentError`.
- `co2lab_compatible`: mesh without connections across pinched-out layers,
  as MRST (see Notes). Default `false`.

# Returns

The NamedTuple `(mesh, columns, grid)`, all new:

- `mesh`: GeoEnergyIO's `UnstructuredMesh`.
- `columns`: its `VEColumns`.
- `grid`: the cut grid section (`Dict` with `"cartDims"`, `"COORD"`,
  `"ZCORN"`, `"ACTNUM"`); its layer 1 is the model's first Johansen layer.

Throws an `ArgumentError` naming the file and its download address if the
deck is not in [`johansen_data_directory`](@ref).

# Notes

- The full-field model has inactive zero-thickness layers between active
  ones. GeoEnergyIO connects cells across them (pinch-out processing, by
  default even without a `PINCH` keyword), so a column continues below.
  MRST does not, and co2lab's `topSurfaceGrid` then drops the cells below
  the gap; `co2lab_compatible = true` gives MRST's behaviour.
"""
function johansen_mesh(model::Symbol = :sector; co2lab_compatible::Bool = false)
    haskey(_JOHANSEN_MODELS, model) ||
        throw(ArgumentError("model must be :sector or :full_field, got :$model"))
    m = _JOHANSEN_MODELS[model]
    path = joinpath(johansen_data_directory(), m.file)
    isfile(path) || throw(ArgumentError("Johansen deck $path not found: download " *
                                        "$(m.url) and unpack it there"))
    # GeoEnergyIO prints a note when it takes the dimensions from the file
    grid = redirect_stdout(devnull) do
        JutulDarcy.GeoEnergyIO.parse_grdecl_file(path)
    end
    grid = _cut_layers(grid, m.layers)
    # co2lab_compatible: no connections across pinched-out layers, as MRST
    mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(
        grid; process_pinch = !co2lab_compatible)
    # cells below holes discarded, as co2lab's topSurfaceGrid
    columns, _ = ve_columns(mesh; below_holes = :discard)
    return (mesh = mesh, columns = columns, grid = grid)
end

# The layers `layers` of a grid section as a new section (co2lab's cutGrdecl):
# ZCORN and ACTNUM cut, COORD kept. Cutting, rather than deactivating the other
# layers, keeps the top faces whole: under inactive cells GeoEnergyIO splits
# them where the layer above is faulted.
function _cut_layers(grid, layers)
    nx, ny, nz = grid["cartDims"]
    # two ZCORN levels per layer (its top and bottom)
    zcorn = reshape(grid["ZCORN"], 2nx, 2ny, 2nz)[:, :, 2first(layers)-1:2last(layers)]
    active = reshape(Bool.(grid["ACTNUM"]), nx, ny, nz)[:, :, layers]
    return Dict{String, Any}("cartDims" => (nx, ny, length(layers)),
                             "COORD" => grid["COORD"], "ZCORN" => vec(zcorn),
                             "ACTNUM" => active)
end

"""
    johansen_well_point() -> NTuple{2, Float64}

The map position of the injection well of the Johansen benchmark
(Eigestad et al.), as co2lab's `trapsJohansen` takes it: the centroid of the
sector model's cell (48, 48) in the formation's top layer.

Returns `(x, y)` (m). Reads the sector model ([`johansen_mesh`](@ref)).
"""
function johansen_well_point()
    m = johansen_mesh(:sector)
    c = findfirst(c -> cell_ijk(m.mesh, c) == (48, 48, 1), 1:number_of_cells(m.mesh))
    x = tpfv_geometry(m.mesh).cell_centroids[:, c]
    return (x[1], x[2])
end

# ---------------------------- The VE sector model -----------------------------

"Files of the sector model read by `johansen_grid_section`."
const _JOHANSEN_VE_FILES = ("NPD5.grdecl", "NPD5_Permeability.txt", "NPD5_Porosity.txt")

"Address of the Johansen data set (MRST's `dataset_johansen`), all files."
const _JOHANSEN_ZIP = "https://www.sintef.no/contentassets/" *
                      "124f261f170947a6bc51dd76aea66129/johansen.zip"

"""
    johansen_grid_section(; min_perm = 0.11) -> Dict{String, Any}

The Johansen sector model (NPD5) with its rock as a GeoEnergyIO grid
section, as co2lab's `makeJohansenVEgrid` prepares it from the data set.
The rock comes from `NPD5_Permeability.txt` (mD) and `NPD5_Porosity.txt`,
one value per cell in deck order.

# Arguments

- `min_perm`: cells with permeability below this (mD) are made inactive;
  the default 0.11 removes the shales (Dunlin above, Amundsen below).

# Returns

A new `Dict{String, Any}` of the deck layers that hold active cells (cut as
co2lab's `cutGrdecl`; deck layers 6–10 for the default `min_perm`):

- `"cartDims"`, `"COORD"`, `"ZCORN"`, `"ACTNUM"`: the grid.
- `"layers"`: the deck layers kept (`UnitRange`); section layer 1 is the
  first of them.
- `"PERMX"`, `"PERMY"` (= x), `"PERMZ"` (0.1 x): permeability (m², as MRST
  converts it), `Array{Float64, 3}`.
- `"PORO"`: porosity, `Array{Float64, 3}`.

Throws an `ArgumentError` naming the directory and the download address if
a file is not in [`johansen_data_directory`](@ref).
"""
function johansen_grid_section(; min_perm::Real = 0.11)
    dir = johansen_data_directory()
    missing_files = filter(f -> !isfile(joinpath(dir, f)), collect(_JOHANSEN_VE_FILES))
    isempty(missing_files) ||
        throw(ArgumentError("Johansen data not found in $dir (missing: " *
                            join(missing_files, ", ") * "): download " *
                            "$_JOHANSEN_ZIP and unpack it there"))
    # GeoEnergyIO prints a note when it takes the dimensions from the file
    g = redirect_stdout(devnull) do
        JutulDarcy.GeoEnergyIO.parse_grdecl_file(joinpath(dir, "NPD5.grdecl"))
    end
    dims = Tuple(g["cartDims"])
    K = _read_johansen_values(joinpath(dir, "NPD5_Permeability.txt"), dims)
    poro = _read_johansen_values(joinpath(dir, "NPD5_Porosity.txt"), dims)

    # the layers holding active cells cut out: under the inactive shales
    # GeoEnergyIO would split the top faces where the layer above is faulted
    active = Bool.(g["ACTNUM"]) .& (K .>= min_perm)
    has_active = [any(active[:, :, k]) for k in axes(active, 3)]
    layers = findfirst(has_active):findlast(has_active)
    grid = _cut_layers(merge(g, Dict("ACTNUM" => active)), layers)

    # m² in MRST's order of operations: [1 1 0.1] .* K .* milli .* darcy
    si(f) = (f .* K[:, :, layers]) .* 1e-3 .* _mrst_darcy()
    return merge!(grid, Dict("PERMX" => si(1.0), "PERMY" => si(1.0),
                             "PERMZ" => si(0.1), "PORO" => poro[:, :, layers],
                             "layers" => layers))
end

# The numbers of a Johansen rock file, in file order (as MRST loads them), as
# a new array of size `dims`; '%' starts a comment
function _read_johansen_values(path, dims)
    v = Float64[]
    for line in eachline(path)
        startswith(lstrip(line), "%") && continue
        append!(v, parse.(Float64, split(line)))
    end
    length(v) == prod(dims) ||
        error("$path has $(length(v)) values, the grid $(prod(dims)) cells")
    return reshape(v, dims)
end

"""
    johansen_ve_domain(; min_perm = 0.11, kwargs...) -> NamedTuple

The Johansen sector model as a VE domain, as co2lab's `exampleVEJohansen`
builds it (`makeJohansenVEgrid` with `remakeGrid`, then `topSurfaceGrid` and
`averageRock`).

Steps: [`johansen_grid_section`](@ref); its largest connected piece (MRST's
`G(1)`), meshed by GeoEnergyIO; columns, with the cells below a hole
discarded; then [`ve_domain`](@ref) with the column averages of x- and
y-permeability and porosity.

# Arguments

- `min_perm`: as [`johansen_grid_section`](@ref) (mD). Default 0.11.
- `kwargs`: passed on to [`ve_domain`](@ref) (e.g. `fault_flux`).

# Returns

The NamedTuple `(mesh, columns, ve, grid, sides)`, all new:

- `mesh`: the 3D `UnstructuredMesh`, discarded cells removed.
- `columns`: its `VEColumns`.
- `ve`: [`ve_domain`](@ref)'s result.
- `grid`: the grid section, `"ACTNUM"` the largest piece (before
  discarding).
- `sides`: the lateral boundary sides ([`ve_boundary_sides`](@ref)), all
  closed until opened by [`ve_boundary_conditions`](@ref).
"""
function johansen_ve_domain(; min_perm::Real = 0.11, kwargs...)
    grid = johansen_grid_section(; min_perm)
    # the largest connected piece only, as co2lab's G(1)
    grid["ACTNUM"] = _largest_piece(grid)
    mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
    # cells below holes dropped, as co2lab's topSurfaceGrid; then the columns
    # of the trimmed mesh
    columns, discarded = ve_columns(mesh; below_holes = :discard)
    mesh = remove_cells(mesh, discarded)
    columns, _ = ve_columns(mesh)

    at(key) = [grid[key][cell_ijk(mesh, c)...] for c in 1:number_of_cells(mesh)]
    ve = ve_domain(mesh, columns; perm_x = at("PERMX"), perm_y = at("PERMY"),
                   poro = at("PORO"), kwargs...)
    sides, _ = ve_boundary_sides(mesh, columns, ve.domain)
    return (mesh = mesh, columns = columns, ve = ve, grid = grid, sides = sides)
end
