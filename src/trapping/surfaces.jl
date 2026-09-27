# Surfaces from files for trap analysis (roadmap Stage 8 item 6): co2lab-
# spillpoint's readIrapClassicAsciiSurf and readIGEMSIRAP — IRAP classic ASCII
# surfaces, and the IGEMS data set's surfaces as one-layer grids — and
# one-layer grids under any surface given on a lattice (synthetic caprocks).
#
# co2lab's IGEMS grids swap the file's x and y (the depth matrix is built with
# y along the grid's x); here the surface keeps the file's orientation, and
# `co2lab_compatible = true` swaps as co2lab does, for comparisons.

# ------------------------------- IRAP surfaces --------------------------------

"""
    read_irap_surface(path) -> NamedTuple

An IRAP classic ASCII surface (co2lab's `readIrapClassicAsciiSurf`): header
of `ny, dx, dy`, `xmin, xmax, ymin, ymax`, `nx, angle`, then `nx · ny` values,
x fastest; the missing value 9999900 becomes `NaN`. Throws if the header's
spacing does not match its extent.

Returns `(x, y, z, angle)`, new: node positions along x and y, `z` an
`nx × ny` matrix (x along the first index — co2lab's is the transpose), and
the rotation angle in degrees (not applied, as in co2lab).
"""
function read_irap_surface(path::AbstractString)
    tok = split(read(path, String))
    val(k) = parse(Float64, tok[k])
    ny, dx, dy = Int(val(2)), val(3), val(4)
    xmin, xmax, ymin, ymax = val(5), val(6), val(7), val(8)
    nx, angle = Int(val(9)), val(10)
    # co2lab checks the spacing exactly
    (xmax - xmin) / (nx - 1) == dx && (ymax - ymin) / (ny - 1) == dy ||
        throw(ArgumentError("$path: inconsistent IRAP header (dx, dy against extent)"))

    data = parse.(Float64, tok[20:end])
    length(data) == nx * ny ||
        throw(ArgumentError("$path: expected $(nx * ny) values, found $(length(data))"))
    z = reshape(replace(data, 9999900.0 => NaN), nx, ny)
    return (x = _linspace(xmin, xmax, nx), y = _linspace(ymin, ymax, ny), z = z,
            angle = angle)
end

# -------------------------------- IGEMS grids ---------------------------------

"""
    igems_surface_path(name, i) -> String

The IGEMS surface file `i` of family `name` (e.g. `"OSSNP1"`, 1–100) in
`atlas_data_directory()`: `IGEMS surfaces/surfaces/<name>/<i>_<name>.irap`,
or for `i == 1` without the full data set the sample
`IGEMSsample/one_of_each/1_<name>.irap` (as co2lab's `readIGEMSIRAP`).
Throws an `ArgumentError` naming the file if it is not there.
"""
function igems_surface_path(name::AbstractString, i::Integer)
    base = atlas_data_directory()
    full = joinpath(base, "IGEMS surfaces", "surfaces", name, "$(i)_$name.irap")
    sample = joinpath(base, "IGEMSsample", "one_of_each", "$(i)_$name.irap")
    isfile(full) && return full
    i == 1 && isfile(sample) && return sample
    throw(ArgumentError("IGEMS surface $full not found: download the IGEMS " *
                        "surfaces (files.nr.no/igems/surfaces.zip) into $base"))
end

"""
    igems_mesh(name, i; coarsening = 1, co2lab_compatible = false) -> NamedTuple

IGEMS surface `i` of family `name` as a one-layer grid (co2lab's
`readIGEMSIRAP(name, i, 'coarse', [c c])`): every `coarsening`-th value
in each direction, then the first and last x-lines replaced by their
neighbours (co2lab's fill of the missing edge values), cells of the
surface's spacing (times `coarsening`) and 100 m thick from the surface
down, cells with an undefined corner removed. With `co2lab_compatible`
the surface is transposed first, as co2lab builds it.

Returns `(mesh, columns)`: a new `UnstructuredMesh` and its `VEColumns`.
"""
function igems_mesh(name::AbstractString, i::Integer; coarsening::Integer = 1,
                    co2lab_compatible::Bool = false)
    s = read_irap_surface(igems_surface_path(name, i))
    z = s.z[1:coarsening:end, 1:coarsening:end]
    # then, as co2lab, the first and last x-lines take their neighbours' values
    z[1, :], z[end, :] = z[2, :], z[end-1, :]
    co2lab_compatible && (z = permutedims(z))
    h = (s.x[2] - s.x[1]) * coarsening
    x, y = [h * i for i in 0:size(z, 1)-1], [h * j for j in 0:size(z, 2)-1]
    return surface_mesh(x, y, z; bottom = z .+ 100.0)
end

# ------------------------------ Surface grids ---------------------------------

"""
    surface_mesh(x, y, top; bottom = top .+ 100.0) -> NamedTuple

A one-layer corner-point grid between two surfaces given at the nodes of a
lattice — for synthetic caprocks (co2lab's examples perturb a `cartGrid`'s
nodes) and read surfaces: `x` (`nx + 1` values) and `y` (`ny + 1`) the node
positions, `top` and `bottom` `(nx + 1) × (ny + 1)` depths (positive down).
Cells with an undefined (`NaN`) corner on `top` are left out.

Returns `(mesh, columns)`: a new `UnstructuredMesh` (GeoEnergyIO's) and its
`VEColumns`.
"""
function surface_mesh(x::AbstractVector, y::AbstractVector, top::AbstractMatrix;
                      bottom::AbstractMatrix = top .+ 100.0)
    nx, ny = length(x) - 1, length(y) - 1
    size(top) == size(bottom) == (nx + 1, ny + 1) ||
        throw(ArgumentError("top and bottom must be $(nx + 1) × $(ny + 1)"))
    defined = isfinite.(top)
    # undefined nodes: a pillar from 0 to 100 m (touching inactive cells only)
    zt, zb = ifelse.(defined, top, 0.0), ifelse.(defined, bottom, 100.0)
    coord = Float64[]
    for j in 0:ny, i in 0:nx
        append!(coord, [x[i+1], y[j+1], zt[i+1, j+1], x[i+1], y[j+1], zb[i+1, j+1]])
    end

    # corner a of 2n along x lies on node 1 + a ÷ 2; layer c: top, bottom
    zcorn = [(c == 1 ? zt : zb)[1 + a ÷ 2, 1 + b ÷ 2]
             for c in 1:2 for b in 1:2ny for a in 1:2nx]
    active = [all(defined[i:i+1, j:j+1]) for i in 1:nx, j in 1:ny, _ in 1:1]
    grid = Dict("cartDims" => (nx, ny, 1), "COORD" => coord, "ZCORN" => zcorn,
                "ACTNUM" => active)
    mesh = JutulDarcy.GeoEnergyIO.mesh_from_grid_section(grid)
    return (mesh = mesh, columns = first(ve_columns(mesh)))
end
