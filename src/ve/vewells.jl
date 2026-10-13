# VE wells (roadmap Stage 2 item 4).
#
# A VE well is one of JutulDarcy's simple wells in one column. As in co2lab
# (`convertwellsVE.m`), it communicates with the full column thickness: its
# well index is the 2D Peaceman index of the column's lateral extent times
# the column height H, and it is supplied directly to `setup_well`.
#
# The pressure datum of the well (the depth its bottom-hole pressure refers to)
# is an option:
# - `:caprock` (co2lab): the well's pressure is compared directly with the
#   column's caprock pressure; no hydrostatic term (co2lab's `dZ = 0`);
# - `:bottom`: the bottom-hole pressure refers to the bottom of the column.
#   JutulDarcy's well model then adds the wellbore hydrostatic term between the
#   bottom and the caprock, where the column's pressure unknown lives.
# The perforation is placed at the caprock in both cases, so the reservoir
# side of the connection is the caprock pressure, as in co2lab.

# --------------------------------- Well index ---------------------------------

# MRST's `wellConstant` table for 'ip_simple' (computeWellIndex.m): aspect ratio
# of the cell against the constant in the equivalent radius
const _WELL_CONSTANT_RATIO = (1, 2, 3, 4, 5, 8, 9, 16, 17, 32, 33, 64, 65)
const _WELL_CONSTANT_VALUE = (0.292, 0.278, 0.262, 0.252, 0.244, 0.231, 0.229,
                              0.220, 0.219, 0.213, 0.213, 0.210, 0.210)

"""
    _mrst_well_constant(d1, d2, inner_product) -> Float64

The constant in the Peaceman equivalent radius, as MRST's `wellConstant`:
0.14 for `:ip_tpf`; for `:ip_simple` a table of the cell's aspect ratio
(rounded as MATLAB's `round`), interpolated and extrapolated linearly.
"""
function _mrst_well_constant(d1, d2, inner_product::Symbol)
    inner_product == :ip_tpf && return 0.14
    inner_product == :ip_simple ||
        throw(ArgumentError("inner_product must be :ip_simple or :ip_tpf, " *
                            "got :$inner_product"))

    # aspect ratio (≥ 1), rounded as MATLAB's round
    r = max(round(d1 / d2, RoundNearestTiesAway), round(d2 / d1, RoundNearestTiesAway))
    x, y = _WELL_CONSTANT_RATIO, _WELL_CONSTANT_VALUE
    # linear in MRST's table, extrapolated beyond its ends
    i = clamp(searchsortedlast(collect(x), r), 1, length(x) - 1)
    return y[i] + (y[i+1] - y[i]) * (r - x[i]) / (x[i+1] - x[i])
end

"""
    ve_well_index(ve, col; radius = 0.1, skin = 0.0, inner_product = :ip_simple)
        -> Float64

co2lab's well index (m³) for a VE well in column `col`: the Peaceman index of
the column's lateral extent (`ve.geometry.dx`, `dy`) with the column's lateral
permeability, times the column height `H` — `convertwellsVE.m`, which computes
the 2D index on the top-surface grid and multiplies it by `H`.

`ve` is the NamedTuple from [`ve_domain`](@ref). `inner_product` selects the
constant in the equivalent radius as MRST does: `:ip_simple` (co2lab's default
in `convertwellsVE`; a table of the cell's aspect ratio) or `:ip_tpf` (0.14,
JutulDarcy's Peaceman constant). The formula itself is JutulDarcy's
`compute_peaceman_index`, which is the same as MRST's. `ve` is only read.
"""
function ve_well_index(ve, col::Integer; radius = 0.1, skin = 0.0,
                       inner_product::Symbol = :ip_simple)
    g = ve.geometry
    dx, dy, H = g.dx[col], g.dy[col], g.H[col]
    k = ve.domain[:permeability][col]
    c = _mrst_well_constant(dx, dy, inner_product)

    # Kh = H·k: co2lab's 2D index (Kh = k) times H
    return JutulDarcy.compute_peaceman_index((dx, dy, H), k, radius;
                                             constant = c, skin = skin)
end

# ----------------------------- Wells and controls -----------------------------

"""
    setup_ve_well(ve, col; name, datum = :caprock, WI = nothing, radius = 0.1,
                  skin = 0.0, inner_product = :ip_simple, kwarg...) -> DataDomain

A VE well in column `col`: a JutulDarcy simple well with one perforation,
for use in [`setup_co2lab_model`](@ref)'s `wells`. `ve` is the NamedTuple from
[`ve_domain`](@ref); it is only read.

- `WI`: the well index (m³). Defaults to co2lab's,
  [`ve_well_index`](@ref)`(ve, col; radius, skin, inner_product)`.
- `datum`: the depth the well's bottom-hole pressure refers to.
  `:caprock` is co2lab's convention: no hydrostatic term between the well and
  the column's caprock pressure (in a fine cell of a hybrid grid: its
  centroid, where its pressure is; the perforation is there too). `:bottom` puts it at the bottom of the column
  (caprock depth + `H`); the well model then adds the wellbore hydrostatic
  term `ρ·g·H` (JutulDarcy's `PerforationGravityDifference` = `-g·H`, with
  the wellbore mixture density). Under rate control the two give the same
  flow; under pressure control the BHP targets differ by that term.
- `kwarg...` go to JutulDarcy's `setup_well` (e.g. `dir`).

Returns a new well `DataDomain` (JutulDarcy's `setup_well`): the perforation
at the column's caprock depth (a fine cell's centroid), the well node at the datum depth, `WI` as its
`:well_index`.
"""
function setup_ve_well(ve, col::Integer; name::Symbol, datum::Symbol = :caprock,
                       WI = nothing, radius = 0.1, skin = 0.0,
                       inner_product::Symbol = :ip_simple, kwarg...)
    WI = something(WI, ve_well_index(ve, col; radius, skin, inner_product))
    return setup_ve_well(ve, [col]; name, datum, WI = [WI], radius, skin, kwarg...)
end

"""
    setup_ve_well(ve, cols::AbstractVector; name, datum = :caprock, WI = nothing,
                  radius = 0.1, skin = 0.0, inner_product = :ip_simple, kwarg...)
        -> DataDomain

A VE well perforating several blocks `cols` (e.g. the stacked columns of one
pillar, `ve_columns(...; below_holes = :stack)`): a JutulDarcy simple well
with one perforation per block, each at the depth of that block's pressure
unknown (a VE column's caprock, a fine cell's centroid). `WI` is one index
per block (default [`ve_well_index`](@ref) of each); `datum = :caprock` puts
the well's pressure datum at the shallowest perforation, `:bottom` at the
bottom of the deepest block. The well model adds the wellbore's hydrostatic
term between the datum and each deeper perforation. Other arguments as for
one column. Returns a new well `DataDomain`; `ve` is only read.
"""
function setup_ve_well(ve, cols::AbstractVector{<:Integer}; name::Symbol,
                       datum::Symbol = :caprock, WI = nothing, radius = 0.1,
                       skin = 0.0, inner_product::Symbol = :ip_simple, kwarg...)
    cols = collect(Int, cols)
    WI = something(WI, [ve_well_index(ve, c; radius, skin, inner_product)
                        for c in cols])
    length(WI) == length(cols) || throw(ArgumentError("one WI per block"))
    z0 = ve.domain[:reference_depth][cols]

    # datum: the shallowest perforation, or the bottom of the deepest block
    g = ve.geometry
    deep = cols[argmax(g.z_top[cols] .+ g.H[cols])]
    z_ref = datum == :caprock ? minimum(z0) : _well_datum_depth(g, deep, datum)

    # perforations at the depth of each block's pressure unknown: a VE
    # column's caprock, a fine cell's centroid
    centers = copy(ve.domain[:cell_centroids])
    centers[3, cols] = z0

    return setup_well(ve.domain, cols; name = name, WI = collect(Float64, WI),
                      radius = radius, skin = skin, cell_centers = centers,
                      reference_depth = z_ref, kwarg...)
end

# depth of a VE well's pressure datum in column `col` (z0: the column's
# pressure datum, its caprock or a fine cell's centroid)
function _well_datum_depth(geometry, col, datum::Symbol; z0 = geometry.z_top[col])
    datum == :caprock && return z0
    datum == :bottom && return geometry.z_top[col] + geometry.H[col]
    throw(ArgumentError("datum must be :caprock or :bottom, got :$datum"))
end

"""
    co2_injector_control(target, fluid::CO2LabVEFluid) -> InjectorControl

JutulDarcy's `InjectorControl` for pure CO2 injection with a co2lab fluid
([`co2lab_ve_fluid`](@ref)): injected mixture `[0, 1]`, "surface" density =
co2lab's CO2 reference density (surface rates are volumes at that density,
not at 1 atm), and the wellbore column of the well's hydrostatic term filled
with CO2 (`phases = ((2, 1.0),)`; JutulDarcy's default is phase 1, brine,
which gives a wrong BHP whenever the datum is not at the perforation).
`fluid` is only read.

Returns a new `InjectorControl` with the given `target` (e.g.
`TotalRateTarget(q)` with `q` in m³/s at the reference density, or
`BottomHolePressureTarget(p)`).
"""
co2_injector_control(target, fluid::CO2LabVEFluid) =
    InjectorControl(target, [0.0, 1.0], density = fluid.rho_ref[2],
                    phases = ((2, 1.0),))
