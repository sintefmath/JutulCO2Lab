# co2lab's VE constitutive functions and PVT (roadmap Stage 2 items 1 and 3).
#
# - the vertical model (roadmap Stage 6 item 0): how a column's fine-scale
#   profile follows from its upscaled state. It gives the heights (h, h_max)
#   from sG, sGmax, the column height H and Δρ, the relative permeabilities
#   from the heights, and the fine-scale saturations for reconstruction. h is
#   the capillary height Pc/(Δρ g): the plume thickness for a sharp interface.
#   `SharpInterface` is co2lab's 'sharp_interface_simple': upscaledSat2height
#   (homogeneous case), and the relperm with endpoints and the accretion layer
#   of addVERelpermSharpInterface — from the caprock rugosity, a property of
#   the VE grid (`ve_domain`), per column;
# - pcWG = Δρ g h for every vertical model, added to the brine
#   (reference-phase) pressure to give the CO2 pressure; Δρ from our black-oil
#   phase densities (researcher's decision: natural in JutulDarcy; equals
#   co2lab's without dissolved CO2);
# - co2lab's PVT (makeVEFluid): densities linear in pressure per phase,
#   constant or linear viscosities.
#
# makeVEFluid also carries two properties that are not the fluid's; here they
# live elsewhere (researcher, 2026-09-25): rugosity with the VE grid
# (`ve_domain(...; rugosity)`), rock compressibility with the model
# (`rock_compressibility.jl`).
#
# The formulas are plain functions (`_plume_heights` etc.), called by the Jutul
# variables and checked directly against co2lab evaluated under Octave
# (scripts/mrst_dump_ve_fluid.m).

# ------------------------------- Vertical models -------------------------------

"""
    VerticalModel

How the fine-scale profile of a VE column follows from its upscaled state
(roadmap Stage 6 item 0). A model given to [`co2lab_ve_fluid`](@ref) is a
grid-independent specification, from which [`setup_co2lab_model`](@ref) builds
the column-specific object (`_vertical_model`), or already built for the grid
([`SharpInterfaceIntegrated`](@ref)). A model implements
`_column_heights`, `_column_relperm`, `_fine_saturations` and
`_category_saturations`; the capillary pressure is `Δρ g h` for all.
"""
abstract type VerticalModel end

"""
    SharpInterface()

co2lab's `'sharp_interface_simple'`: a sharp interface between the plume and
the brine below it, in a vertically homogeneous column, with residual CO2
between `h` and `h_max` and the accretion layer of the caprock rugosity.
"""
struct SharpInterface <: VerticalModel end

"""
    SharpInterfaceIntegrated

co2lab's `'sharp_interface_integrated'`: a sharp interface in a column whose
rock varies with depth. Heights by pore volume, relperm by permeability, both
integrated over the column's fine cells; built with
[`SharpInterfaceIntegrated(mesh, columns; poro, perm_x)`](@ref
SharpInterfaceIntegrated(::Any, ::VEColumns)). Per column, at the boundaries
of its fine cells (top first, `pos[c]:pos[c+1]-1`), the depth below the
caprock and the integrals from the caprock down.
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

The column-specific vertical model for the reservoir `domain` (a
`DataDomain` from `ve_domain`) and `fluid`, built from the specification
`spec`. Models without column data return `spec` itself. Arguments are only
read.
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

co2lab's VE fluid, with `makeVEFluid`'s option names, so an MRST case
translates one to one:
- `vertical`: the [`VerticalModel`](@ref), makeVEFluid's relperm model —
  [`SharpInterface`](@ref) for `'sharp_interface_simple'`;
- `residual = (res_water, res_gas)`; `krmax = (krw, krg)`: `krw` the brine
  relperm in the residual-CO2 zone (brine below it has relperm 1), `krg` the
  CO2 relperm in the plume; a negative value means one minus the other
  phase's residual saturation;
- `*_rho_pvt`, `*_mu_pvt`: each phase's density and viscosity —
  `:table` (default): JutulDarcy's CO2–brine tables in pressure and
  temperature (see [`jutuldarcy_property_tables`](@ref), configured by
  `tables`, e.g. `tables = (co2_source = :salo24, salt_names = ["NaCl"],
  salt_mole_fractions = [0.01])`); a [`PhaseTable`](@ref) (e.g. co2lab's,
  [`co2lab_sampled_table`](@ref)); `(c, p_ref)` for co2lab's linear
  `ref·(1 + c·(p − p_ref))`; `nothing` for constant. Tabulated properties
  take the temperature of each column (`ve_domain(...; temperature)`).
  `default_pvt` is what a property not given means: `:table`, or `nothing`
  (co2lab's constant; as before Stage 7, for ports of co2lab cases).
- `*_rho_ref`: reference densities, the black-oil "surface" densities
  (`b = ρ/ρ_ref`). For a linear or constant density the density at `p_ref`
  (default co2lab's: CO2 760, brine 1100 kg/m³). For a tabulated one only the
  scale of `b` and of the dissolved-CO2 ratio (below): default JutulDarcy's
  surface densities (CO2 1.868, brine 998.2 kg/m³), or the table's densities
  at `reference_conditions = (p, T)`.
- `*_mu_ref`: constant or linear viscosities only (default co2lab's: CO2
  6e-5, brine 8e-4 Pa·s).
- `dissolution`, `dis_rate`, `dis_max`: CO2 dissolution into brine, with
  `dis_max` the saturated dissolved-CO2 ratio `rs`: the volume of the
  dissolved CO2 at `co2_rho_ref` per volume of brine at `wat_rho_ref`
  (black-oil "surface" volumes, co2lab's convention), i.e. a mass ratio of
  `rs · co2_rho_ref / wat_rho_ref` kg CO2 per kg brine. Or give that mass
  ratio as `dis_max_mass`. Default: co2lab's, `dis_max = 0.03` at its
  reference densities, ≈ 0.0207 kg/kg (whatever the reference densities).
  `dis_rate = 0` or `Inf`: instantaneous (brine saturated wherever CO2 is
  present); otherwise rate-driven, with co2lab's transfer rule
  (`computeDissolutionTransfer`).

With `default_pvt = :table`, a reference value given without its `*_pvt` is
refused: state `*_pvt` (`nothing` for co2lab's constant property). Returns a
new `CO2LabVEFluid`, with the endpoints, properties and reference values
resolved.
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
                         dissolution, Float64(dis_rate), dmax)
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

Settings of co2lab's VE fluid; build with [`co2lab_ve_fluid`](@ref). Phase
order in all tuples: (brine, CO2). Type parameters: the vertical model, and
the density and viscosity specifications. `rho_pvt`/`mu_pvt` hold per phase
`(c, p_ref)` (linear; `(0.0, 0.0)` for constant) or a [`PhaseTable`](@ref).
`dissolution`, `dis_rate`, `dis_max` as makeVEFluid's options.
"""
struct CO2LabVEFluid{M <: VerticalModel, R, U}
    vertical::M                  # vertical model (specification)
    res_water::Float64           # residual brine saturation
    res_gas::Float64             # residual CO2 saturation
    krw::Float64                 # brine relperm in the residual-CO2 zone
    krg::Float64                 # CO2 relperm in the plume
    rho_ref::NTuple{2, Float64}  # reference ("surface") densities
    rho_pvt::R                   # densities: (c, p_ref) or PhaseTable, per phase
    mu_ref::NTuple{2, Float64}   # reference viscosities (NaN for a table)
    mu_pvt::U                    # viscosities, as rho_pvt
    dissolution::Bool            # CO2 dissolves into brine
    dis_rate::Float64            # dissolution rate; 0 or Inf: instantaneous
    dis_max::Float64             # saturated rs
end

"""
    instantaneous_dissolution(fluid) -> Bool

Whether `fluid` has dissolution with co2lab's instantaneous mode (`dis_rate`
0 or `Inf`).
"""
instantaneous_dissolution(f::CO2LabVEFluid) =
    f.dissolution && (f.dis_rate == 0 || isinf(f.dis_rate))

# `nothing` → constant; (c, p_ref) → linear; makeVEFluid's 4-entry table spec
# is refused (its tables: co2lab_sampled_table).
function _linear_pvt(spec)
    isnothing(spec) && return (0.0, 0.0)
    length(spec) == 2 ||
        error("makeVEFluid's table range ($spec) is not taken; give :table " *
              "(JutulDarcy's), co2lab_sampled_table(file), nothing (constant) " *
              "or (c, p_ref) (linear)")
    return (Float64(spec[1]), Float64(spec[2]))
end

# ---------------------------------- Formulas ----------------------------------

# Clamps with MRST's tie rule (ADI.m, max/min): against a constant, a tie keeps
# the variable and its derivative; between two variables, a tie keeps the
# second. Julia's max(x, 0) and min(u, v) break ties the other way, which only
# changes derivatives at exact ties — but a cell that has just switched to
# saturated sits at sG = 0 exactly, and there a lost derivative makes Newton
# overshoot (instantaneous dissolution cycled without these).
"`max(x, c)` for a constant `c`; a tie keeps `x`."
_max_c(x, c) = x < c ? zero(x) + c : x
"`min(x, c)` for a constant `c`; a tie keeps `x`."
_min_c(x, c) = x > c ? zero(x) + c : x
"`min(u, v)` of two variables; a tie keeps `v`."
_min_v(u, v) = u < v ? u : v

"""
    _plume_heights(sG, sGmax, H, rw, rg) -> (h, h_max)

co2lab's `upscaledSat2height` for a vertically homogeneous column
(`upscaledSat2height.m`, `sharp_interface_case`): `h_max = sGmax·H/(1−rw)`,
`h = (H·sG − h_max·rg)/(1−rg−rw)`, clamped to `[0, H]`.
"""
function _plume_heights(sG, sGmax, H, rw, rg)
    h_max = sGmax * H / (1 - rw)
    h = (H * sG - h_max * rg) / (1 - rg - rw)
    return (_min_c(_max_c(h, 0.0), H), h_max)
end

"""
    _sharp_interface_relperm(h, h_max, H, krw, krg, dh) -> (krW, krG)

`addVERelpermSharpInterface.m`, `krW_simple`/`krG_simple`: CO2 flows in the
plume above the accretion layer `krH_lost = min(dh, H)`; brine in the
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
model `m`, from the maximum CO2 saturation `smax` (as it stands within the
step), the maximum capillary height `h_max`, the column height `H` and the
density difference `drho` (co2lab's `sW_untouched`, for the dissolution
transfer). Default: the column below `h_max`, `(H − min(h_max, H))/H`.
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

Brine minus CO2 density at pressure `p`, dissolved-CO2 ratio `rs` and
temperature `T`, as the model's black-oil densities (`update_bo_densities!`):
the brine carries its dissolved CO2.
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
([`ve_domain`](@ref)); zero where the domain has none. The relperm loses the
accretion layer `min(dh, H)` under the caprock.
"""
struct Rugosity <: Jutul.ScalarVariable end

function Jutul.default_parameter_values(data_domain, model, ::Rugosity, symb)
    haskey(data_domain, :rugosity) && return copy(data_domain[:rugosity])
    return zeros(Jutul.number_of_cells(Jutul.physical_representation(data_domain)))
end

# ------------------------ VE / fine-cell branch points ------------------------
#
# The three VE variables below branch per cell on `is_ve` (`nothing` = all VE;
# else one Bool per cell), so that roadmap Stage 12 (hybrid: fine 3D cells)
# adds the fine side without reshaping them
# (analysis/ve-meshing/isve-branch-proposal.md). With `nothing` the branch is
# compiled away. Saturation regions, when they come, are a separate `regions`
# field in JutulDarcy's form, not part of this flag.

"Is cell `i` a VE column? Always for `is_ve = nothing`."
_is_ve(::Nothing, i) = true
_is_ve(is_ve::AbstractVector{Bool}, i) = is_ve[i]

# Fine-cell sides of the relperm and capillary pressure: roadmap Stage 12.
_fine_relperm(fine, i) = error("fine 3D cells (relperm) are roadmap Stage 12")
_fine_capillary_pressure(i) = error("fine 3D cells (capillary pressure) are " *
                                    "roadmap Stage 12")

"""
Capillary heights per column (m): row 1 `h`, row 2 `h_max` — the plume
thickness for a sharp interface, `Pc/(Δρ g)` in general — from the vertical
model `vertical`. Zero in fine cells (`is_ve`, see above), where nothing
reads them.
"""
struct PlumeHeights{M <: VerticalModel, V} <: Jutul.VectorVariables
    vertical::M         # vertical model (column-specific)
    res_water::Float64  # residual brine saturation
    res_gas::Float64    # residual CO2 saturation
    is_ve::V            # true for VE columns, per cell
end
PlumeHeights(vertical, rw, rg; is_ve = nothing) =
    PlumeHeights(vertical, Float64(rw), Float64(rg), is_ve)
Jutul.values_per_entity(model, ::PlumeHeights) = 2

# Uses sGmax_eff = max(min(sGmax, 1 − rw), sG): the stored sGmax is only
# updated after a converged step, so within a step the current sG must count
# too (MRST's updateState does this after every Newton update,
# CO2VEBlackOilTypeModel.m:352-353). Otherwise h_max < h while the plume grows.
# `drho(i)`: the density difference of cell i, or nothing.
# Modifies: `hh` (plume heights).
function _update_heights!(hh, v::PlumeHeights, sG, sGmax, H, drho, ix)
    rw, rg = v.res_water, v.res_gas
    for i in ix
        if !_is_ve(v.is_ve, i)
            hh[1, i], hh[2, i] = zero(eltype(hh)), zero(eltype(hh))
            continue
        end

        # sGmax as it stands within this step (see above)
        smax = max(min(sGmax[i], 1 - rw), sG[i])
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
Upscaled relative permeabilities (brine, CO2) of the VE columns from the
vertical model `vertical`. In fine cells (`is_ve`) the 3D relperm `fine` —
roadmap Stage 12, `nothing` until then.
"""
struct VERelPerm{M <: VerticalModel, V, F} <: JutulDarcy.AbstractRelativePermeabilities
    vertical::M         # vertical model (column-specific)
    krw::Float64        # brine relperm in the residual-CO2 zone
    krg::Float64        # CO2 relperm in the plume
    res_water::Float64  # residual brine saturation
    res_gas::Float64    # residual CO2 saturation
    is_ve::V            # true for VE columns, per cell
    fine::F             # 3D relperm of fine cells
end
VERelPerm(vertical, krw, krg, rw, rg; is_ve = nothing, fine = nothing) =
    VERelPerm(vertical, Float64.((krw, krg, rw, rg))..., is_ve, fine)

# `sat(i)`: (sG, sGmax_eff) of cell i, or nothing; `drho(i)`: the density
# difference of cell i, or nothing.
# Modifies: `kr` (relative permeabilities).
function _update_relperm!(kr, v::VERelPerm, hh, H, dh, sat, drho, ix)
    for i in ix
        if !_is_ve(v.is_ve, i)
            kr[1, i], kr[2, i] = _fine_relperm(v.fine, i)
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
    sat(i) = (sG[i], max(min(sGmax[i], 1 - rw), sG[i]))
    _update_relperm!(kr, v, PlumeHeights, ColumnHeight, Rugosity, sat,
                     i -> ρ[1, i] - ρ[2, i], ix)
end

# The sharp interfaces: heights only, no saturations or densities.
# Modifies: `kr` (relative permeabilities).
@jutul_secondary function update_ve_relperm!(kr, v::VERelPerm{<:_DensityFreeModel},
                                            model, PlumeHeights, ColumnHeight, Rugosity,
                                            ix)
    _update_relperm!(kr, v, PlumeHeights, ColumnHeight, Rugosity, i -> nothing,
                     i -> nothing, ix)
end

"""
`pcWG = (ρ_brine − ρ_CO2)·g·h`, added to the brine pressure for CO2. Fine cells
(`is_ve`): roadmap Stage 12.
"""
struct VECapillaryPressure{V} <: JutulDarcy.AbstractCapillaryPressure
    is_ve::V  # true for VE columns, per cell
end
VECapillaryPressure(; is_ve = nothing) = VECapillaryPressure(is_ve)

# Modifies: `pc` (capillary pressure, one row: CO2 relative to brine).
@jutul_secondary function update_ve_pc!(pc, v::VECapillaryPressure, model, PlumeHeights,
                                       PhaseMassDensities, ix)
    g = Jutul.gravity_constant
    for i in ix
        if !_is_ve(v.is_ve, i)
            pc[1, i] = _fine_capillary_pressure(i)
            continue
        end

        # the brine–CO2 density difference over the plume height
        pc[1, i] = (PhaseMassDensities[1, i] - PhaseMassDensities[2, i]) * g *
                   PlumeHeights[1, i]
    end
end

"Viscosities `μ_ref·(1 + c·(p − p_ref))` per phase (brine, CO2)."
struct LinearViscosities <: JutulDarcy.PhaseVariables
    mu_ref::NTuple{2, Float64}          # reference viscosities
    pvt::NTuple{2, NTuple{2, Float64}}  # (c, p_ref) per phase
end

# Modifies: `mu` (phase viscosities).
@jutul_secondary function update_linear_visc!(mu, v::LinearViscosities, model, Pressure,
                                             ix)
    for i in ix, ph in 1:2
        mu[ph, i] = _linear(v.mu_ref[ph], v.pvt[ph], Pressure[i])
    end
end

# --------------------- Installing the VE fluid in a model ---------------------

"""
    _set_ve_fluid!(model, fluid, m = _vertical_model(fluid.vertical, ...))

Install co2lab's VE constitutive functions in the reservoir submodel `model`
(which must be on a domain from `ve_domain`), with the column-specific vertical
model `m` (built from `fluid` unless given).

Modifies: `model.parameters` (adds `:ColumnHeight`, `:Rugosity`),
`model.secondary_variables` (`:PlumeHeights`, `:RelativePermeabilities`,
`:CapillaryPressure`), `model.output_variables` (adds `:PlumeHeights`,
`:CapillaryPressure`). `fluid` and `m` are only read.
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

    # the column-specific vertical model is shared by heights and relperm
    rw, rg = fluid.res_water, fluid.res_gas
    S[:PlumeHeights] = PlumeHeights(m, rw, rg; is_ve)
    S[:RelativePermeabilities] = VERelPerm(m, fluid.krw, fluid.krg, rw, rg; is_ve)
    S[:CapillaryPressure] = VECapillaryPressure(; is_ve)

    push!(model.output_variables, :PlumeHeights, :CapillaryPressure)
    unique!(model.output_variables)
    return model
end

"""
    _co2lab_pvt(fluid) -> (shrinkage, viscosity)

The fluid's PVT as Jutul variables: with tabulated properties `PVTShrinkage`
and `PVTViscosities` (pressure and temperature); otherwise co2lab's linear
`b = 1 + c·(p − p_ref)` per phase (the reference density is the "surface"
density) and linear or constant viscosities (pressure only). Returns new
variable objects; `fluid` is only read.
"""
function _co2lab_pvt(fluid::CO2LabVEFluid)
    rho, mu = fluid.rho_pvt, fluid.mu_pvt
    if _tabulated(rho...)
        b = PVTShrinkage(fluid.rho_ref, rho)
    else
        (cw, pw), (cg, pg) = rho
        b = LinearShrinkage((1.0, 1.0), (cw, cg), (pw, pg))
    end
    visc = _tabulated(mu...) ? PVTViscosities(fluid.mu_ref, mu) :
                               LinearViscosities(fluid.mu_ref, mu)
    return (b, visc)
end

"""
    needs_temperature(fluid) -> Bool

Whether `fluid` has a tabulated property, and so needs each column's
temperature (`ve_domain(...; temperature)`; 30 °C where none is given).
"""
needs_temperature(f::CO2LabVEFluid) = _tabulated(f.rho_pvt..., f.mu_pvt...)
