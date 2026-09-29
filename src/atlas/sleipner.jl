# The Sleipner benchmark model (IEAGHG, layer 9 of the Utsira formation at
# Sleipner) as a VE domain (roadmap Stage 5 item 4) — co2lab-common/
# makeSleipnerVEmodel.m, as calibrateSleipner calls it ('assign_coords', true).
#
# The data set is no longer downloadable; it is read from a local copy in
# `sleipner_data_directory()`. Unlike co2lab's default, which negates y and so
# mirrors the grid in local coordinates, the grid is placed in its global UTM
# coordinates (zone 31, like the Atlas) from the file's MAPAXES.

# ---------------------------------- The data ----------------------------------

"The files of the IEAGHG Sleipner benchmark (layer 9) that are read."
const _SLEIPNER_FILES = ("M9X1.grdecl", "M9X1_perm_X_mD_.inc", "M9X1_perm_Z_mD_.inc",
                         "M9X1_poro___.inc")

"""
    sleipner_data_directory() -> String

Directory holding the Sleipner benchmark files (`M9X1.grdecl` and the
`M9X1_*.inc` rock files): `Sleipner` in [`atlas_data_directory`](@ref), as
MRST's data directory has it.
"""
sleipner_data_directory() = joinpath(atlas_data_directory(), "Sleipner")

# Throw an informative error if the benchmark files are not there
function _check_sleipner_data(dir)
    missing_files = filter(f -> !isfile(joinpath(dir, f)), collect(_SLEIPNER_FILES))
    isempty(missing_files) && return nothing
    throw(ArgumentError(
        "Sleipner benchmark data not found in $dir (missing: " *
        join(missing_files, ", ") * "). The IEAGHG Sleipner benchmark model is no " *
        "longer publicly downloadable: place a copy of its files there (or set " *
        "CO2LAB_ATLAS_DATA to the directory that contains `Sleipner`)."))
end

# ------------------------------ The grid section ------------------------------

"""
    sleipner_grid_section(; min_perm = 200.0) -> Dict{String, Any}

The Sleipner benchmark's grid as a GeoEnergyIO grid section, as
co2lab's `makeSleipnerVEmodel('assign_coords', true)` prepares it: the grid
of `M9X1.grdecl` moved to global coordinates by its `MAPAXES` (co2lab's
`mapAxes`), cells with x-permeability below `min_perm` (mD) made inactive
(the shale cells), and the rock of the `M9X1_*.inc` files — `"PERMX"`,
`"PERMY"` (= x) and `"PERMZ"` in m² (MRST's `convertFrom(perm, milli*darcy)`)
and `"PORO"`, as `Array{Float64, 3}`. Throws an `ArgumentError` naming the
directory if the data are not present.

Returns a new `Dict` with `"cartDims"`, `"COORD"`, `"ZCORN"`, `"ACTNUM"` and
the rock arrays.
"""
function sleipner_grid_section(; min_perm::Real = 200.0)
    dir = sleipner_data_directory()
    _check_sleipner_data(dir)
    # GeoEnergyIO prints a note when it takes the dimensions from the file
    g = redirect_stdout(devnull) do
        JutulDarcy.GeoEnergyIO.parse_grdecl_file(joinpath(dir, "M9X1.grdecl"))
    end
    dims = Tuple(g["cartDims"])
    rock(file, key) = JutulDarcy.GeoEnergyIO.parse_grdecl_file(joinpath(dir, file), dims;
                                                               input_units = :si)[key]
    permx = rock("M9X1_perm_X_mD_.inc", "PERMX")
    permz = rock("M9X1_perm_Z_mD_.inc", "PERMZ")

    md = 1e-3 * _mrst_darcy()                          # milli*darcy, MRST's order
    # shale cells (x-permeability below min_perm) inactive; permeability in m²
    return Dict{String, Any}("cartDims" => dims,
        "COORD" => _map_axes(g["COORD"], g["MAPAXES"]), "ZCORN" => g["ZCORN"],
        "ACTNUM" => g["ACTNUM"] .& (permx .>= min_perm),
        "PERMX" => permx .* md, "PERMY" => permx .* md, "PERMZ" => permz .* md,
        "PORO" => rock("M9X1_poro___.inc", "PORO"))
end

# COORD (x, y, z per pillar end) moved to global coordinates: MRST's mapAxes,
# x e1 + y e2 + origin, with the axes and origin of MAPAXES (y-axis point,
# origin, x-axis point)
function _map_axes(coord, mapaxes)
    o = mapaxes[3:4]
    e1 = (mapaxes[5:6] .- o) ./ hypot((mapaxes[5:6] .- o)...)
    e2 = (mapaxes[1:2] .- o) ./ hypot((mapaxes[1:2] .- o)...)
    c = reshape(copy(coord), 3, :)
    for p in axes(c, 2)
        x, y = c[1, p], c[2, p]
        c[1, p] = (x * e1[1] + y * e2[1]) + o[1]
        c[2, p] = (x * e1[2] + y * e2[2]) + o[2]
    end
    return vec(c)
end

# -------------------------------- The VE domain -------------------------------

"""
    sleipner_ve_domain(; min_perm = 200.0) -> NamedTuple

The Sleipner benchmark as a VE domain (co2lab's `makeSleipnerVEmodel` with
global coordinates): [`sleipner_grid_section`](@ref), GeoEnergyIO's mesh,
columns with cells below a hole discarded (co2lab's `topSurfaceGrid`
default), and [`ve_domain`](@ref) with the column averages of x-permeability
and porosity (`averageRock`). All lateral sides of the domain are open in
co2lab's Sleipner examples (`bcIxVE`): `sides` are the domain's
[`ve_boundary_sides`](@ref).

Returns the NamedTuple `(mesh, columns, ve, grid, sides)`, all new: the 3D
`UnstructuredMesh` (discarded cells removed), its `VEColumns`,
[`ve_domain`](@ref)'s result, the grid section, and the boundary sides.
"""
function sleipner_ve_domain(; min_perm::Real = 200.0)
    grid = sleipner_grid_section(; min_perm)
    mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
    # cells below holes dropped, then the columns of the trimmed mesh
    columns, discarded = ve_columns(mesh)
    mesh = ve_trim(mesh, discarded)
    columns, _ = ve_columns(mesh)

    at(key) = [grid[key][cell_ijk(mesh, c)...] for c in 1:number_of_cells(mesh)]
    ve = ve_domain(mesh, columns; perm_x = at("PERMX"), poro = at("PORO"))
    sides, _ = ve_boundary_sides(mesh, columns, ve.domain)
    return (mesh = mesh, columns = columns, ve = ve, grid = grid, sides = sides)
end
