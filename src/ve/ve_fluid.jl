# co2lab's VE fluid: settings, constitutive functions and their Jutul variables.
#
# - Vertical models: how a column's fine-scale profile follows from its
#   upscaled state.  They give the capillary heights (h, h_max) from sG, sGmax,
#   the column height H and Δρ; the relperms from the heights; and the
#   fine-scale saturations for reconstruction.  h = Pc/(Δρ g), the plume
#   thickness for a sharp interface.
# - Capillary pressure pcWG = Δρ g h for every vertical model, added to the
#   brine (reference-phase) pressure to give the CO2 pressure.  Δρ comes from
#   the black-oil phase densities (equal to co2lab's without dissolved CO2).
# - PVT: co2lab's (`makeVEFluid`: densities linear in pressure, constant or
#   linear viscosities) or JutulDarcy's CO2–brine tables.
#
# Two of makeVEFluid's options are not the fluid's and live elsewhere: the
# rugosity with the VE grid (`ve_domain(...; rugosity)`), rock compressibility
# with the model (`setup_co2lab_model(...; rock_compressibility)`).
#
# The formulas are plain functions (`_plume_heights` etc.), called by the Jutul
# variables and testable directly against co2lab.

# ------------------------------- Vertical models -------------------------------

"""
    VerticalModel

The family of VE vertical models: how the fine-scale profile of a column
follows from its upscaled state (saturations, column height, Δρ).  It is
given to [`co2lab_ve_fluid`](@ref) as `vertical`.

A model is either a grid-independent specification, from which
[`setup_co2lab_model`](@ref) builds the column-specific object
(`_vertical_model`), or already built for a grid (e.g.
[`SharpInterfaceIntegrated`](@ref), [`LayeredFringe`](@ref)).

# Subtypes

- [`SharpInterface`](@ref): homogeneous column, sharp interface.
- [`SharpInterfaceIntegrated`](@ref): sharp interface, rock varying with depth.
- [`CapillaryFringeTable`](@ref): capillary fringe in a homogeneous column.
- [`LayeredFringe`](@ref): capillary fringe, rock varying with depth.

# Notes

- A model implements `_column_heights`, `_column_relperm`,
  `_fine_saturations` and `_category_saturations`.
- The capillary pressure is `Δρ g h` for all models.
"""
abstract type VerticalModel end

"""
    SharpInterface()

co2lab's `'sharp_interface_simple'` vertical model (the default): a sharp
interface between the plume and the brine below it, in a vertically
homogeneous column.

- Heights as co2lab's `upscaledSat2height` (homogeneous case); residual CO2
  lies between `h` and `h_max`.
- Relperm as co2lab's `addVERelpermSharpInterface`: endpoints `krw`, `krg`,
  and no CO2 flow in the accretion layer `min(dh, H)` under the caprock (CO2
  entering it is trapped; the plume flows below), with `dh` the rugosity per
  column ([`ve_domain`](@ref)).
"""
struct SharpInterface <: VerticalModel end

"""
    SharpInterfaceIntegrated

co2lab's `'sharp_interface_integrated'` vertical model: a sharp interface in
a column whose rock varies with depth.  Heights follow from pore volume,
relperm from permeability, both integrated over the column's fine cells.

Built for one grid with
[`SharpInterfaceIntegrated(mesh, columns; poro, perm_x)`](@ref
SharpInterfaceIntegrated(::Any, ::VEColumns)).

# Fields

Per column `c`, one entry per boundary of its fine cells, top first, at
`pos[c]:pos[c+1]-1`: the depth below the caprock and the integrals from the
caprock down to that boundary.
"""
struct SharpInterfaceIntegrated <: VerticalModel
    pos::Vector{Int}    # column c: entries pos[c]:pos[c+1]-1, one per cell boundary
    z::Vector{Float64}  # depth of the boundary below the caprock (m)
    V::Vector{Float64}  # pore volume above it, ∫ φ·area dz (m³)
    K::Vector{Float64}  # x-permeability integrated down to it, ∫ k_x dz (m³)
end

# Vertical models whose heights and relperm depend on the saturations and the
# column alone (not on the densities)
const _DensityFreeModel = Union{SharpInterface, SharpInterfaceIntegrated}

"""
    _vertical_model(spec, domain, fluid) -> VerticalModel

Build the column-specific vertical model from the specification `spec`, for
the reservoir `domain` (a `DataDomain` from [`ve_domain`](@ref)) and `fluid`.
Returns `spec` itself for models without column data.
"""
_vertical_model(spec::VerticalModel, domain, fluid) = spec

# ---------------------------------- Settings ----------------------------------

"""
    co2lab_ve_fluid(; vertical = SharpInterface(), residual = (0.0, 0.0),
                      krmax = (-1.0, -1.0),
                      co2_rho_pvt = :table, wat_rho_pvt = :table,
                      co2_mu_pvt = :table, wat_mu_pvt = :table,
                      co2_rho_ref = nothing, wat_rho_ref = nothing,
                      co2_mu_ref = nothing, wat_mu_ref = nothing,
                      default_pvt = :table, tables = (;),
                      reference_conditions = nothing, dissolution = false,
                      dis_rate = 0.0, dis_max = nothing, dis_max_mass = nothing)
        -> CO2LabVEFluid

Build co2lab's VE fluid: the vertical model, relperm endpoints, PVT and
dissolution.  The keywords are `makeVEFluid`'s options, so an MRST case
translates one to one.  The result goes to [`setup_co2lab_model`](@ref).

# Arguments

- `vertical`: the [`VerticalModel`](@ref) (makeVEFluid's relperm model);
  default [`SharpInterface`](@ref), co2lab's `'sharp_interface_simple'`.
- `residual = (res_water, res_gas)`: residual brine and CO2 saturations.
- `krmax = (krw, krg)`: relperm endpoints; `krw` the brine relperm in the
  residual-CO2 zone (brine below it has relperm 1), `krg` the CO2 relperm in
  the plume.  A negative value means one minus the other phase's residual
  saturation.
- `co2_rho_pvt`, `wat_rho_pvt`, `co2_mu_pvt`, `wat_mu_pvt`: each phase's
  density and viscosity, one of
  - `:table`: JutulDarcy's CO2–brine tables in pressure and temperature
    ([`jutuldarcy_property_tables`](@ref)), configured by `tables`;
  - a [`PhaseTable`](@ref) (e.g. co2lab's sampled tables, read from an MRST export);
  - `(c, p_ref)`: co2lab's linear `ref·(1 + c·(p − p_ref))`;
  - `nothing`: constant `ref`.
- `default_pvt`: what a property not given means: `:table` (default), or
  `nothing` (co2lab's constant, for ports of co2lab cases).
- `tables`: keywords for [`jutuldarcy_property_tables`](@ref), e.g.
  `(co2_source = :salo24, salt_names = ["NaCl"], salt_mole_fractions = [0.01])`.
- `co2_rho_ref`, `wat_rho_ref`: reference densities (kg/m³), the black-oil
  "surface" densities (`b = ρ/ρ_ref`).
  - Linear or constant density: the density at `p_ref`; default co2lab's,
    CO2 760, brine 1100.
  - Tabulated density: only the scale of `b` and of `rs`; default
    JutulDarcy's surface densities (CO2 1.868, brine 998.2), or the table's
    densities at `reference_conditions = (p, T)` if given.
- `co2_mu_ref`, `wat_mu_ref`: reference viscosities (Pa·s), for constant or
  linear viscosities only; default co2lab's, CO2 6e-5, brine 8e-4.
- `dissolution`: whether CO2 dissolves into brine.
- `dis_rate`: `0` or `Inf` for instantaneous dissolution (brine saturated
  wherever CO2 is present); otherwise rate-driven, with co2lab's transfer
  rule (`computeDissolutionTransfer`).
- `dis_max`: saturated dissolved-CO2 ratio `rs`, in black-oil "surface"
  volumes (co2lab's convention): volume of dissolved CO2 at `co2_rho_ref` per
  volume of brine at `wat_rho_ref`, i.e. `rs · co2_rho_ref / wat_rho_ref` kg
  CO2 per kg brine.
- `dis_max_mass`: alternatively, that mass ratio (kg/kg); not with `dis_max`.
  Default for both: co2lab's `dis_max = 0.03` at its reference densities,
  ≈ 0.0207 kg/kg whatever the reference densities.

Tabulated properties take each column's temperature
(`ve_domain(...; temperature)`).

# Returns

A new [`CO2LabVEFluid`](@ref), with endpoints, properties and reference
values resolved.

# Notes

- With `default_pvt = :table`, a reference value given without its `*_pvt`
  is an error: state `*_pvt` (`nothing` for co2lab's constant property).
- makeVEFluid's rugosity goes to [`ve_domain`](@ref), its rock
  compressibility (`pvMult_*`) to [`setup_co2lab_model`](@ref).
"""
function co2lab_ve_fluid(; vertical::VerticalModel = SharpInterface(),
                           residual = (0.0, 0.0), krmax = (-1.0, -1.0),
                           co2_rho_pvt = missing, wat_rho_pvt = missing,
                           co2_mu_pvt = missing, wat_mu_pvt = missing,
                           co2_rho_ref = nothing, wat_rho_ref = nothing,
                           co2_mu_ref = nothing, wat_mu_ref = nothing,
                           default_pvt = :table, tables = (;),
                           reference_conditions = nothing,
                           dissolution::Bool = false, dis_rate = 0.0, dis_max = nothing,
                           dis_max_mass = nothing)
    # makeVEFluid's rugosity goes to ve_domain, pvMult_* to setup_co2lab_model
    rw, rg = Float64.(residual)
    krw = krmax[1] < 0 ? 1 - rg : Float64(krmax[1])
    krg = krmax[2] < 0 ? 1 - rw : Float64(krmax[2])

    # each phase property: tabulated (JutulDarcy's by default) or co2lab's
    jd = Ref{Any}(nothing)
    jd_tables() = isnothing(jd[]) ? (jd[] = jutuldarcy_property_tables(; tables...)) :
                  jd[]
    spec(x, ref, name, jd) = _pvt_spec(x, ref, name, jd, default_pvt)
    rho = (spec(wat_rho_pvt, wat_rho_ref, "wat_rho", () -> jd_tables().density[1]),
           spec(co2_rho_pvt, co2_rho_ref, "co2_rho", () -> jd_tables().density[2]))
    mu = (spec(wat_mu_pvt, wat_mu_ref, "wat_mu", () -> jd_tables().viscosity[1]),
          spec(co2_mu_pvt, co2_mu_ref, "co2_mu", () -> jd_tables().viscosity[2]))

    rho_ref = _reference_densities(rho, (wat_rho_ref, co2_rho_ref), reference_conditions)
    mu_ref = _reference_viscosities(mu, (wat_mu_ref, co2_mu_ref))
    dmax = _dis_max(dis_max, dis_max_mass, rho_ref)
    return CO2LabVEFluid(vertical, rw, rg, krw, krg, rho_ref, rho, mu_ref, mu,
                         dissolution, Float64(dis_rate), dmax, nothing)
end

# One phase property: `missing` (not given) is `default`, a table only if its
# reference value is not given; `:table` JutulDarcy's (`jd()`); `nothing` or
# `(c, p_ref)` linear
function _pvt_spec(spec, ref, name, jd, default)
    if ismissing(spec) && default !== :table
        spec = default
    elseif ismissing(spec)
        isnothing(ref) ||
            error("$(name)_ref given without $(name)_pvt: state $(name)_pvt " *
                  "(nothing for constant, (c, p_ref), :table or a PhaseTable)")
        spec = :table
    end
    spec === :table && return jd()
    spec isa PhaseTable && return spec
    return _linear_pvt(spec)
end

# Reference densities: given; for a table its densities at
# `conditions = (p, T)`, else JutulDarcy's surface densities; for a linear or
# constant density co2lab's defaults
function _reference_densities(rho, given, conditions)
    surface, co2lab = (998.207150, 1.868048), (1100.0, 760.0)
    return ntuple(2) do ph
        isnothing(given[ph]) || return Float64(given[ph])
        rho[ph] isa PhaseTable || return co2lab[ph]
        isnothing(conditions) ? surface[ph] : Float64(rho[ph](conditions...))
    end
end

# Reference viscosities: co2lab's defaults; none (NaN) for a table
function _reference_viscosities(mu, given)
    return ntuple(2) do ph
        if mu[ph] isa PhaseTable
            isnothing(given[ph]) || error("a reference viscosity is given for a " *
                                          "tabulated viscosity")
            return NaN
        end
        return isnothing(given[ph]) ? (8e-4, 6e-5)[ph] : Float64(given[ph])
    end
end

# The saturated dissolved-CO2 ratio rs: given, from a mass ratio, or co2lab's
# default capacity (0.03 at its reference densities 760, 1100)
function _dis_max(dis_max, dis_max_mass, (rhoW, rhoG))
    isnothing(dis_max) || isnothing(dis_max_mass) ||
        error("give dis_max or dis_max_mass, not both")
    isnothing(dis_max) || return Float64(dis_max)
    mass = something(dis_max_mass, 0.03 * 760.0 / 1100.0)
    # a mass ratio m_CO2/m_brine is rs·ρ_CO2,ref/ρ_brine,ref
    (rhoW, rhoG) == (1100.0, 760.0) && isnothing(dis_max_mass) && return 0.03
    return Float64(mass) * rhoW / rhoG
end

"""
    CO2LabVEFluid

Settings of co2lab's VE fluid, with all defaults resolved; build with
[`co2lab_ve_fluid`](@ref).  Phase order in all tuples: (brine, CO2).

# Fields

- `rho_pvt`, `mu_pvt`: per phase `(c, p_ref)` (linear; `(0.0, 0.0)` for
  constant) or a [`PhaseTable`](@ref).
- `dissolution`, `dis_rate`, `dis_max`: as the keywords of
  [`co2lab_ve_fluid`](@ref).
- `fine`: the fluid of a hybrid grid's fine 3D cells (a `FineScaleFluid` with
  its rock type per cell, from [`co2lab_ve_fluid`](@ref)`(::FineScaleFluid,
  …)`), or `nothing`.

Type parameters: the vertical model, the density and the viscosity
specifications, and the fine cells' fluid.
"""
struct CO2LabVEFluid{M <: VerticalModel, R, U, F}
    vertical::M                  # vertical model (specification)
    res_water::Float64           # residual brine saturation
    res_gas::Float64             # residual CO2 saturation
    krw::Float64                 # brine relperm in the residual-CO2 zone
    krg::Float64                 # CO2 relperm in the plume
    rho_ref::NTuple{2, Float64}  # reference ("surface") densities (kg/m³)
    rho_pvt::R                   # densities: (c, p_ref) or PhaseTable, per phase
    mu_ref::NTuple{2, Float64}   # reference viscosities (Pa·s; NaN for a table)
    mu_pvt::U                    # viscosities, as rho_pvt
    dissolution::Bool            # CO2 dissolves into brine
    dis_rate::Float64            # dissolution rate; 0 or Inf: instantaneous
    dis_max::Float64             # saturated dissolved-CO2 ratio rs
    fine::F                      # fine cells' fluid (hybrid grids), or nothing
end

"""
    instantaneous_dissolution(fluid) -> Bool

Returns `true` if `fluid` dissolves CO2 instantaneously (co2lab's mode for
`dis_rate` 0 or `Inf`), `false` otherwise, including without dissolution.
"""
instantaneous_dissolution(f::CO2LabVEFluid) =
    f.dissolution && (f.dis_rate == 0 || isinf(f.dis_rate))

# `nothing` → constant; (c, p_ref) → linear; makeVEFluid's 4-entry table range
# is refused (give its tables as a PhaseTable instead).
function _linear_pvt(spec)
    isnothing(spec) && return (0.0, 0.0)
    length(spec) == 2 ||
        error("makeVEFluid's table range ($spec) is not taken; give :table " *
              "(JutulDarcy's), a PhaseTable, nothing (constant) " *
              "or (c, p_ref) (linear)")
    return (Float64(spec[1]), Float64(spec[2]))
end

# ---------------------------------- Formulas ----------------------------------

# The plain value of `x`, every level of AD derivatives removed: what branches
# decide on.  ForwardDiff (from 1.0) breaks a tie between a value with
# derivatives and another by their derivatives; co2lab's rules decide on the
# values, and a tie must take the branch the rule names.
_plain(x::Jutul.ForwardDiff.Dual) = _plain(Jutul.ForwardDiff.value(x))
_plain(x) = x

# Clamps with MRST's tie rule for AD max/min: against a constant, a tie keeps
# the variable and its derivative; between two variables, the second.  Julia's
# max(x, 0) and min(u, v) break ties the other way.  That matters: a cell just
# switched to saturated sits at sG = 0 exactly, and a lost derivative there
# makes Newton overshoot (instantaneous dissolution cycles).
"`max(x, c)` for a constant `c`; a tie keeps `x`."
_max_c(x, c) = _plain(x) < _plain(c) ? zero(x) + c : x
"`min(x, c)` for a constant `c`; a tie keeps `x`."
_min_c(x, c) = _plain(x) > _plain(c) ? zero(x) + c : x
"`min(u, v)` of two variables; a tie keeps `v`."
_min_v(u, v) = _plain(u) < _plain(v) ? u : v

"""
    _effective_max(sGmax, sG, rw, rg) -> Real

`sGmax` as it stands within a step: `max(min(sGmax, 1 − rw), sG)`.

- With hysteresis (`rg > 0`), `sGmax` is an unknown, projected to at least
  `sG` after every Newton update (as co2lab's).  A tie keeps `sGmax`, so this
  is `sGmax` itself with co2lab's Jacobian; otherwise Newton can cycle where
  a cell switches between drainage and imbibition.
- Without, `sGmax` is a parameter updated after converged steps, so the
  current `sG` must count while the plume grows; a tie keeps `sG` and its
  derivative.
"""
function _effective_max(sGmax, sG, rw, rg)
    m = _min_c(sGmax, 1 - rw)
    rg > 0 && return _plain(m) >= _plain(sG) ? m : sG
    return _plain(m) > _plain(sG) ? m : sG
end

"""
    _plume_heights(sG, sGmax, H, rw, rg) -> (h, h_max)

co2lab's `upscaledSat2height` for a vertically homogeneous sharp-interface
column: `h_max = sGmax·H/(1−rw)`, `h = (H·sG − h_max·rg)/(1−rg−rw)`, with `h`
clamped to `[0, H]`.
"""
function _plume_heights(sG, sGmax, H, rw, rg)
    h_max = sGmax * H / (1 - rw)
    h = (H * sG - h_max * rg) / (1 - rg - rw)
    return (_min_c(_max_c(h, 0.0), H), h_max)
end

"""
    _sharp_interface_relperm(h, h_max, H, krw, krg, dh) -> (krW, krG)

co2lab's `addVERelpermSharpInterface` (`krW_simple`, `krG_simple`): CO2 flows
in the plume below the accretion layer `min(dh, H)`; brine flows in the
brine-only zone below `h_max`, and with endpoint `krw` in the residual-CO2
zone between `h` and `h_max`.
"""
function _sharp_interface_relperm(h, h_max, H, krw, krg, dh)
    lost = min(dh, H)
    # CO2 flows in the plume, below the accretion layer
    krG = _max_c(h - lost, 0.0) * krg / H
    # brine: relperm 1 below h_max, krw in the residual-CO2 zone above it
    krW = (H - _max_c(h_max, lost) + _max_c(h_max - _max_c(h, lost), 0.0) * krw) / H
    return (krW, krG)
end

# ------------------------- Vertical-model interface ---------------------------

"""
    _column_heights(m, i, sG, sGmax, H, drho, rw, rg) -> (h, h_max)

Capillary heights of column `i` under vertical model `m`, from its CO2
saturation `sG`, maximum CO2 saturation `sGmax`, height `H`, brine–CO2
density difference `drho` and residual saturations `rw`, `rg`.
"""
_column_heights(::SharpInterface, i, sG, sGmax, H, drho, rw, rg) =
    _plume_heights(sG, sGmax, H, rw, rg)

"""
    _column_relperm(m, i, sG, sGmax, h, h_max, H, drho, dh, krw, krg, rw, rg)
        -> (krW, krG)

Upscaled relative permeabilities (brine, CO2) of column `i` under vertical
model `m`, from its saturations, heights, column height `H`, density
difference `drho`, rugosity `dh`, endpoints `krw`, `krg` and residual
saturations `rw`, `rg`.
"""
_column_relperm(::SharpInterface, i, sG, sGmax, h, h_max, H, drho, dh, krw, krg, rw,
                rg) = _sharp_interface_relperm(h, h_max, H, krw, krg, dh)

"""
    _fine_saturations(m, i, zb, h, h_max, drho, rw, rg) -> (s, smax)

CO2 saturation and maximum CO2 saturation of each fine cell of column `i`
(cumulative cell-bottom depths below the caprock `zb`) under vertical model
`m`, with density difference `drho`. Returns two new vectors, top cell
first.
"""
_fine_saturations(::SharpInterface, i, zb, h, h_max, drho, rw, rg) =
    _column_saturations(zb, h, h_max, rw, rg)

"""
    _untouched_brine(m, i, smax, h_max, H, drho) -> Real

Fraction of column `i`'s pore volume that CO2 has never reached under vertical
model `m` (co2lab's `sW_untouched`, for the dissolution transfer), from the
maximum CO2 saturation `smax` as it stands within the step, the maximum
capillary height `h_max`, the column height `H` and the density difference
`drho`.  Default: the column below `h_max`, `(H − min(h_max, H))/H`.
"""
_untouched_brine(::VerticalModel, i, smax, h_max, H, drho) = (H - _min_c(h_max, H)) / H

"Linear property `ref·(1 + c·(p − p_ref))`; constant when `c = 0`."
_linear(ref, (c, p_ref), p) = ref * (1 + c * (p - p_ref))

"""
    _shrinkage(fluid, p, T) -> (bW, bG)

The fluid's shrinkage factors `b = ρ(p, T)/ρ_ref` (brine, CO2) at pressure `p`
and temperature `T` (ignored by linear properties).
"""
_shrinkage(fluid::CO2LabVEFluid, p, T) =
    ntuple(ph -> _phase_value(fluid.rho_ref[ph], fluid.rho_pvt[ph], p, T) /
                 fluid.rho_ref[ph], 2)

"""
    _density_difference(fluid, p, rs, T = nothing) -> Real

Brine minus CO2 density (kg/m³) at pressure `p`, dissolved-CO2 ratio `rs` and
temperature `T`, consistent with the model's black-oil densities: the brine
carries its dissolved CO2.
"""
function _density_difference(fluid::CO2LabVEFluid, p, rs, T = nothing)
    (ρW, ρG), (bW, bG) = fluid.rho_ref, _shrinkage(fluid, p, T)
    return bW * (ρW + rs * ρG) - bG * ρG
end

# ------------------------------ Jutul variables -------------------------------

"Column height `H` per VE column (m), read from `domain[:column_height]`."
struct ColumnHeight <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, ::ColumnHeight, symb)
    haskey(data_domain, :column_height) ||
        error("the domain has no :column_height; build it with ve_domain")
    return copy(data_domain[:column_height])
end

"""
Caprock rugosity `dh` per VE column (m), read from `domain[:rugosity]`
([`ve_domain`](@ref)); zero if the domain has none.  The relperm loses the
accretion layer `min(dh, H)` under the caprock: CO2 entering it is trapped.
"""
struct Rugosity <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, ::Rugosity, symb)
    haskey(data_domain, :rugosity) && return copy(data_domain[:rugosity])
    return zeros(Jutul.number_of_cells(Jutul.physical_representation(data_domain)))
end

# ------------------------ VE / fine-cell branch points ------------------------
#
# The three VE variables below branch per cell on `is_ve`: `nothing` when all
# cells are VE columns (the branch is then compiled away), else one Bool per
# cell, false for the fine 3D cells of a hybrid grid.  Saturation regions, if
# added, belong in a separate `regions` field in JutulDarcy's form, not here.

"Is cell `i` a VE column? Always for `is_ve = nothing`."
_is_ve(::Nothing, i) = true
_is_ve(is_ve::AbstractVector{Bool}, i) = is_ve[i]

# The fine cells' relperm and capillary pressure come from the fluid's
# `FineCellFluid` (`fine_fluid.jl`).

"""
Capillary heights per column (m), from the vertical model: row 1 `h`, row 2
`h_max`.  For a sharp interface the plume thickness, `Pc/(Δρ g)` in general.
Zero in fine cells (`is_ve` false), where nothing reads them.
"""
struct PlumeHeights{M <: VerticalModel, V} <: Jutul.VectorVariables
    vertical::M         # vertical model (column-specific)
    res_water::Float64  # residual brine saturation
    res_gas::Float64    # residual CO2 saturation
    is_ve::V            # per cell, true for VE columns; nothing: all VE
end
PlumeHeights(vertical, rw, rg; is_ve = nothing) =
    PlumeHeights(vertical, Float64(rw), Float64(rg), is_ve)
Jutul.values_per_entity(model, ::PlumeHeights) = 2

# Heights from sGmax as it stands within the step (`_effective_max`): the
# stored sGmax is updated only after a converged step, so the current sG must
# count too (as co2lab does after every Newton update); otherwise h_max < h
# while the plume grows.  `drho(i)`: cell i's density difference, or nothing.
# Modifies: `hh` (plume heights).
function _update_heights!(hh, v::PlumeHeights, sG, sGmax, H, drho, ix)
    rw, rg = v.res_water, v.res_gas
    for i in ix
        if !_is_ve(v.is_ve, i)
            hh[1, i], hh[2, i] = zero(eltype(hh)), zero(eltype(hh))
            continue
        end

        # sGmax as it stands within this step (see above)
        smax = _effective_max(sGmax[i], sG[i], rw, rg)
        hh[1, i], hh[2, i] = _column_heights(v.vertical, i, sG[i], smax, H[i], drho(i),
                                             rw, rg)
    end
    return hh
end

# Any vertical model: the heights may depend on the density difference.
# Modifies: `hh` (plume heights).
@jutul_secondary function update_plume_heights!(hh, v::PlumeHeights, model, sG, sGmax,
                                               ColumnHeight, PhaseMassDensities, ix)
    ρ = PhaseMassDensities
    _update_heights!(hh, v, sG, sGmax, ColumnHeight, i -> ρ[1, i] - ρ[2, i], ix)
end

# The sharp interfaces do not read the densities, so they do not depend on
# them (an unchanged dependency graph and derivatives).
# Modifies: `hh` (plume heights).
@jutul_secondary function update_plume_heights!(hh, v::PlumeHeights{<:_DensityFreeModel},
                                               model, sG, sGmax, ColumnHeight, ix)
    _update_heights!(hh, v, sG, sGmax, ColumnHeight, i -> nothing, ix)
end

"""
Upscaled relative permeabilities (brine, CO2) of the VE columns, from the
vertical model.  In fine cells (`is_ve` false), the 3D relperm `fine`.
"""
struct VERelPerm{M <: VerticalModel, V, F} <: JutulDarcy.AbstractRelativePermeabilities
    vertical::M         # vertical model (column-specific)
    krw::Float64        # brine relperm in the residual-CO2 zone
    krg::Float64        # CO2 relperm in the plume
    res_water::Float64  # residual brine saturation
    res_gas::Float64    # residual CO2 saturation
    is_ve::V            # per cell, true for VE columns; nothing: all VE
    fine::F             # 3D relperm of fine cells (FineCellFluid), or nothing
end
VERelPerm(vertical, krw, krg, rw, rg; is_ve = nothing, fine = nothing) =
    VERelPerm(vertical, Float64.((krw, krg, rw, rg))..., is_ve, fine)

# `sat(i)`: (sG, sGmax_eff) of cell i, or nothing; `drho(i)`: the density
# difference of cell i, or nothing.
# Modifies: `kr` (relative permeabilities).
function _update_relperm!(kr, v::VERelPerm, hh, H, dh, sat, drho, ix)
    for i in ix
        if !_is_ve(v.is_ve, i)
            kr[1, i], kr[2, i] = _fine_relperm(v.fine, i, sat(i)...)
            continue
        end

        s, smax = something(sat(i), (nothing, nothing))
        kr[1, i], kr[2, i] = _column_relperm(v.vertical, i, s, smax, hh[1, i], hh[2, i],
                                             H[i], drho(i), dh[i], v.krw, v.krg,
                                             v.res_water, v.res_gas)
    end
    return kr
end

# Any vertical model: the relperm may depend on the saturations and the
# density difference. sGmax as it stands within the step (as for the heights).
# Modifies: `kr` (relative permeabilities).
@jutul_secondary function update_ve_relperm!(kr, v::VERelPerm, model, sG, sGmax,
                                            PlumeHeights, ColumnHeight, Rugosity,
                                            PhaseMassDensities, ix)
    ρ, rw = PhaseMassDensities, v.res_water
    sat(i) = (sG[i], _effective_max(sGmax[i], sG[i], rw, v.res_gas))
    _update_relperm!(kr, v, PlumeHeights, ColumnHeight, Rugosity, sat,
                     i -> ρ[1, i] - ρ[2, i], ix)
end

# The sharp interfaces: heights only, no saturations or densities (without
# fine cells, which read their saturations).
# Modifies: `kr` (relative permeabilities).
@jutul_secondary function update_ve_relperm!(kr,
                                            v::VERelPerm{<:_DensityFreeModel, Nothing},
                                            model, PlumeHeights, ColumnHeight, Rugosity,
                                            ix)
    _update_relperm!(kr, v, PlumeHeights, ColumnHeight, Rugosity, i -> nothing,
                     i -> nothing, ix)
end

"""
VE capillary pressure `pcWG = (ρ_brine − ρ_CO2)·g·h` (Pa), added to the brine
pressure to give the CO2 pressure.  In fine cells (`is_ve` false), the
capillary pressure of `fine`.
"""
struct VECapillaryPressure{V, F} <: JutulDarcy.AbstractCapillaryPressure
    is_ve::V  # per cell, true for VE columns; nothing: all VE
    fine::F   # fine cells' fluid (FineCellFluid), or nothing
end
VECapillaryPressure(; is_ve = nothing, fine = nothing) = VECapillaryPressure(is_ve, fine)

# The VE columns' pcWG: the brine–CO2 density difference over the plume height
_ve_pc(PlumeHeights, ρ, i) = (ρ[1, i] - ρ[2, i]) * Jutul.gravity_constant *
                             PlumeHeights[1, i]

# All VE: heights and densities only.
# Modifies: `pc` (capillary pressure, one row: CO2 relative to brine).
@jutul_secondary function update_ve_pc!(pc, v::VECapillaryPressure{Nothing}, model,
                                       PlumeHeights, PhaseMassDensities, ix)
    for i in ix
        pc[1, i] = _ve_pc(PlumeHeights, PhaseMassDensities, i)
    end
end

# With fine cells, which read their saturations (sGmax as it stands within the
# step, as for the heights).
# Modifies: `pc` (capillary pressure, one row: CO2 relative to brine).
@jutul_secondary function update_ve_pc!(pc, v::VECapillaryPressure, model, PlumeHeights,
                                       PhaseMassDensities, sG, sGmax, ix)
    f = v.fine
    for i in ix
        if _is_ve(v.is_ve, i)
            pc[1, i] = _ve_pc(PlumeHeights, PhaseMassDensities, i)
        else
            smax = _effective_max(sGmax[i], sG[i], f.res_water, f.res_gas)
            pc[1, i] = _fine_capillary_pressure(f, i, sG[i], smax)
        end
    end
end

"""
co2lab's linear (or constant) viscosities `μ_ref·(1 + c·(p − p_ref))` per
phase (brine, CO2).

- `P = true` (reservoir): at the phase pressures, as co2lab's viscosities
  (unlike its densities); CO2 at `p + pcWG`.
- `P = false` (wells, which have no capillary pressure): at `p`.

Built by `LinearViscosities(mu_ref, pvt; phase_pressures = true)`.
"""
struct LinearViscosities{P} <: JutulDarcy.PhaseVariables
    mu_ref::NTuple{2, Float64}          # reference viscosities (Pa·s)
    pvt::NTuple{2, NTuple{2, Float64}}  # (c, p_ref) per phase; c = 0: constant
end
LinearViscosities(mu_ref, pvt; phase_pressures::Bool = true) =
    LinearViscosities{phase_pressures}(mu_ref, pvt)

# Modifies: `mu` (phase viscosities).
@jutul_secondary function update_linear_visc!(mu, v::LinearViscosities{true}, model,
                                             Pressure, CapillaryPressure, ix)
    for i in ix
        mu[1, i] = _linear(v.mu_ref[1], v.pvt[1], Pressure[i])
        mu[2, i] = _linear(v.mu_ref[2], v.pvt[2], Pressure[i] + CapillaryPressure[1, i])
    end
end

# Modifies: `mu` (phase viscosities).
@jutul_secondary function update_linear_visc!(mu, v::LinearViscosities{false}, model,
                                             Pressure, ix)
    for i in ix, ph in 1:2
        mu[ph, i] = _linear(v.mu_ref[ph], v.pvt[ph], Pressure[i])
    end
end

# --------------------- Installing the VE fluid in a model ---------------------

# A fluid's fine-scale part, required on a grid with fine cells
function _fine_spec(fluid)
    isnothing(fluid.fine) &&
        error("fine 3D cells need a fluid built from a fine-scale one: " *
              "co2lab_ve_fluid(::FineScaleFluid, mesh, columns; ...)")
    return fluid.fine
end

"""
    _set_ve_fluid!(model, fluid, m = _vertical_model(fluid.vertical, ...))

Install co2lab's VE constitutive functions in the reservoir submodel `model`,
which must be on a domain from [`ve_domain`](@ref).  `m` is the
column-specific vertical model, built from `fluid` unless given.  Returns
`model` itself.

Modifies: `model.parameters` (adds `:ColumnHeight`, `:Rugosity`),
`model.secondary_variables` (`:PlumeHeights`, `:RelativePermeabilities`,
`:CapillaryPressure`), `model.primary_variables` (with fine cells `:sG`
becomes a [`FineGasSaturation`](@ref)), `model.output_variables` (adds
`:PlumeHeights`, `:CapillaryPressure`).
"""
function _set_ve_fluid!(model, fluid::CO2LabVEFluid,
                       m = _vertical_model(fluid.vertical, model.data_domain, fluid))
    haskey(model.data_domain, :column_height) ||
        error("VE fluid needs domain[:column_height]; build the domain with ve_domain")

    prm, S = model.parameters, model.secondary_variables
    prm[:ColumnHeight] = ColumnHeight()
    prm[:Rugosity] = Rugosity()

    # per-cell VE flag; `nothing` (branch compiled away) when all cells are VE
    d = model.data_domain
    is_ve = haskey(d, :is_ve) && !all(d[:is_ve]) ? collect(Bool, d[:is_ve]) : nothing
    # fine cells: their functions from the fine-scale fluid
    fine = isnothing(is_ve) ? nothing : _fine_cell_fluid(_fine_spec(fluid), d)

    # the column-specific vertical model is shared by heights and relperm
    rw, rg = fluid.res_water, fluid.res_gas
    S[:PlumeHeights] = PlumeHeights(m, rw, rg; is_ve)
    S[:RelativePermeabilities] = VERelPerm(m, fluid.krw, fluid.krg, rw, rg; is_ve, fine)
    S[:CapillaryPressure] = VECapillaryPressure(; is_ve, fine)

    # fine cells: Newton's CO2 saturation update stopped at the curves' kinks
    P = model.primary_variables
    if !isnothing(fine) && get(P, :sG, nothing) isa GasSaturation
        P[:sG] = FineGasSaturation(_saturation_kinks(_fine_spec(fluid).fluid),
                                   fine.rock_type)
    end

    push!(model.output_variables, :PlumeHeights, :CapillaryPressure)
    unique!(model.output_variables)
    return model
end

"""
    _co2lab_pvt(fluid) -> (shrinkage, viscosity, well_viscosity)

The fluid's PVT as new Jutul variables: the shrinkage factors, the
reservoir's viscosities (at the phase pressures, as co2lab's) and the wells'
(at the pressure).

- Tabulated: `PVTShrinkage`, `PVTViscosities` (pressure and temperature).
- Otherwise co2lab's linear `b = 1 + c·(p − p_ref)` per phase (the reference
  density is the "surface" density), and linear or constant viscosities.
"""
function _co2lab_pvt(fluid::CO2LabVEFluid)
    rho, mu = fluid.rho_pvt, fluid.mu_pvt
    if _tabulated(rho...)
        b = PVTShrinkage(fluid.rho_ref, rho)
    else
        (cw, pw), (cg, pg) = rho
        b = LinearShrinkage((1.0, 1.0), (cw, cg), (pw, pg))
    end
    # the reservoir's at the phase pressures (co2lab's), the wells' at p
    visc(pp) = _tabulated(mu...) ?
               PVTViscosities(fluid.mu_ref, mu; phase_pressures = pp) :
               LinearViscosities(fluid.mu_ref, mu; phase_pressures = pp)
    return (b, visc(true), visc(false))
end

"""
    needs_temperature(fluid) -> Bool

Returns `true` if `fluid` has a tabulated property and so needs each column's
temperature (`ve_domain(...; temperature)`; 30 °C if none is given).
"""
needs_temperature(f::CO2LabVEFluid) = _tabulated(f.rho_pvt..., f.mu_pvt...)
