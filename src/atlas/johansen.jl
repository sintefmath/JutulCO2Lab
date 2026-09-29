# The Johansen formation's NPD models — the sector model (NPD5) and the
# full-field model — as co2lab-spillpoint's trapsJohansen example reads them:
# the layers of the Johansen formation cut from the decks. Both are faulted
# corner-point grids.

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
`FULLFIELD_IMAXJMAX.GRDECL` and its fault include file): `Johansen` in
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
