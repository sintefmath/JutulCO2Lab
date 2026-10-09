# The Johansen formation's NPD models — the sector model (NPD5) and the
# full-field model — as co2lab-spillpoint's trapsJohansen example reads them:
# the layers of the Johansen formation cut from the decks. Both are faulted
# corner-point grids. And the sector model as a VE domain with its rock, as
# co2lab-ve's exampleVEJohansen builds it (makeJohansenVEgrid).

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
model's rock (`NPD5_Permeability.txt`, `NPD5_Porosity.txt`): `Johansen` in
[`atlas_data_directory`](@ref), as MRST's data directory has it.
"""
johansen_data_directory() = joinpath(atlas_data_directory(), "Johansen")

# ------------------------------- The grids ------------------------------------

"""
    johansen_mesh(model = :sector; co2lab_compatible = false) -> NamedTuple

The Johansen formation in one of the NPD models, `model = :sector` (NPD5,
layers 6–11) or `:full_field` (layers 10–14), cut from the decks as
co2lab's `trapsJohansen` cuts them (`cutGrdecl`).

The full-field model has inactive layers of zero thickness between active
ones. GeoEnergyIO connects the cells across them (its pinch-out processing,
by default also without a `PINCH` keyword), so a column continues below
them; MRST's `processGRDECL` does not, and co2lab's `topSurfaceGrid` then
drops the cells below the gap. `co2lab_compatible = true` meshes without
pinch-out connections, as MRST. Throws an `ArgumentError` naming the file
and where to download it if the deck is not in
[`johansen_data_directory`](@ref).

Returns `(mesh, columns, grid)`, all new: GeoEnergyIO's `UnstructuredMesh`,
its `VEColumns` and the cut grid section (a `Dict`: `"cartDims"`, `"COORD"`,
`"ZCORN"`, `"ACTNUM"`; layer 1 is the model's first Johansen layer).
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
    return (mesh = mesh, columns = first(ve_columns(mesh)), grid = grid)
end

# The layers `layers` of a grid section as a section of their own (co2lab's
# cutGrdecl): ZCORN and ACTNUM cut, COORD kept. Cutting, rather than making
# the other layers inactive, keeps the top faces whole: under inactive cells
# GeoEnergyIO splits them where the layer above is faulted.
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

The map position `(x, y)` of the injection well of the Johansen benchmark
(Eigestad et al.): the centroid of the sector model's cell (48, 48) in the
formation's top layer, as co2lab's `trapsJohansen` takes it.
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

The Johansen sector model (NPD5, all 11 layers) as a GeoEnergyIO grid
section, as co2lab's `makeJohansenVEgrid` prepares it from the data set: the
grid of `NPD5.grdecl`; the rock of `NPD5_Permeability.txt` (mD) and
`NPD5_Porosity.txt`, one value per cell in deck order; cells with
permeability below `min_perm` (mD) made inactive, which removes the shales
(Dunlin above, Amundsen below). The layers holding active cells are cut out
of the deck (co2lab's `cutGrdecl`; for `min_perm = 0.11`, deck layers 6–10),
which leaves the cells as they are but keeps their top faces whole. Throws
an `ArgumentError` naming the directory and the download address if a file
is not in [`johansen_data_directory`](@ref).

Returns a new `Dict` with `"cartDims"`, `"COORD"`, `"ZCORN"`, `"ACTNUM"`,
`"layers"` (the deck layers kept, a `UnitRange`; layer 1 of the section is
the first of them) and the rock as `Array{Float64, 3}`: `"PERMX"`, `"PERMY"`
(= x) and `"PERMZ"` (0.1 x) in m², as MRST converts them, and `"PORO"`.
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

# The numbers of a Johansen rock file in file order (MRST's
# reshape(load(file)', [], 1)) as an array of size `dims`; '%' starts a comment
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

The Johansen sector model as a VE domain, built from the data set as
co2lab's `exampleVEJohansen` builds it (`makeJohansenVEgrid` with
`remakeGrid`, then `topSurfaceGrid` and `averageRock`):
[`johansen_grid_section`](@ref); its largest connected piece (MRST's `G(1)`)
meshed by GeoEnergyIO; columns with the cells below a hole discarded
(`discard_below_holes`); and [`ve_domain`](@ref) with the column averages of
the x- and y-permeability and porosity. `kwargs` are passed on to
[`ve_domain`](@ref) (e.g. `fault_flux`).

Returns the NamedTuple `(mesh, columns, ve, grid, sides)`, all new: the 3D
`UnstructuredMesh` (discarded cells removed), its `VEColumns`,
[`ve_domain`](@ref)'s result, the grid section (`ACTNUM` the largest piece,
before discarding), and the lateral boundary sides
([`ve_boundary_sides`](@ref)), all closed until opened by
[`ve_boundary_conditions`](@ref).
"""
function johansen_ve_domain(; min_perm::Real = 0.11, kwargs...)
    grid = johansen_grid_section(; min_perm)
    # the largest connected piece only, as co2lab's G(1)
    grid["ACTNUM"] = _largest_piece(grid)
    mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
    # cells below holes dropped, then the columns of the trimmed mesh
    columns, discarded = ve_columns(mesh)
    mesh = ve_trim(mesh, discarded)
    columns, _ = ve_columns(mesh)

    at(key) = [grid[key][cell_ijk(mesh, c)...] for c in 1:number_of_cells(mesh)]
    ve = ve_domain(mesh, columns; perm_x = at("PERMX"), perm_y = at("PERMY"),
                   poro = at("PORO"), kwargs...)
    sides, _ = ve_boundary_sides(mesh, columns, ve.domain)
    return (mesh = mesh, columns = columns, ve = ve, grid = grid, sides = sides)
end
