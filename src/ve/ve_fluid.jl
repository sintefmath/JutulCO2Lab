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
grid-independent specification; [[`setup_co2lab_model`](@ref) builds the
column-specific object from it (`_vertical_model`). A model implements
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
                      krmax = (-1.0, -1.0), co2_rho_ref = 760.0, wat_rho_ref = 1100.0,
                      co2_rho_pvt = nothing, wat_rho_pvt = nothing,
                      co2_mu_ref = 6e-5, wat_mu_ref = 8e-4,
                      co2_mu_pvt = nothing, wat_mu_pvt = nothing,
                      dissolution = false, dis_rate = 0.0, dis_max = 0.03)
        -> CO2LabVEFluid

co2lab's VE fluid, with `makeVEFluid`'s option names and defaults, so an MRST
case translates one to one:
- `vertical`: the [`VerticalModel`](@ref), makeVEFluid's relperm model —
  [`SharpInterface`](@ref) for `'sharp_interface_simple'`;
- `residual = (res_water, res_gas)`; `krmax = (krw, krg)`: `krw` the brine
  relperm in the residual-CO2 zone (brine below it has relperm 1), `krg` the
  CO2 relperm in the plume; a negative value means one minus the other
  phase's residual saturation;
- `*_rho_ref`: reference densities, which are also the black-oil "surface"
  densities (so `b = ρ(p)/ρ_ref`);
- `*_rho_pvt`, `*_mu_pvt`: `nothing` (constant) or `(c, p_ref)` for
  `ref·(1 + c·(p − p_ref))`. co2lab's third option, a sampled EOS table (and
  its default for CO2 density), is not supported.
- `dissolution`, `dis_rate`, `dis_max`: CO2 dissolution into brine, with
  `dis_max` the saturated dissolved-CO2 ratio `rs` (surface volume of CO2 per
  surface volume of brine). `dis_rate = 0` or `Inf`: instantaneous (brine
  saturated wherever CO2 is present); otherwise rate-driven, with co2lab's
  transfer rule (`computeDissolutionTransfer`).

Returns a new `CO2LabVEFluid`, with the endpoints resolved.
"""
function co2lab_ve_fluid(; vertical::VerticalModel = SharpInterface(),
                           residual = (0.0, 0.0), krmax = (-1.0, -1.0),
                           co2_rho_ref = 760.0, wat_rho_ref = 1100.0,
                           co2_rho_pvt = nothing, wat_rho_pvt = nothing,
                           co2_mu_ref = 6e-5, wat_mu_ref = 8e-4,
                           co2_mu_pvt = nothing, wat_mu_pvt = nothing,
                           dissolution::Bool = false, dis_rate = 0.0, dis_max = 0.03)
    # makeVEFluid's rugosity goes to ve_domain, pvMult_* to setup_co2lab_model
    rw, rg = Float64.(residual)
    krw = krmax[1] < 0 ? 1 - rg : Float64(krmax[1])
    krg = krmax[2] < 0 ? 1 - rw : Float64(krmax[2])
    lin = _linear_pvt

    return CO2LabVEFluid(vertical, rw, rg, krw, krg,
                         (Float64(wat_rho_ref), Float64(co2_rho_ref)),
                         (lin(wat_rho_pvt), lin(co2_rho_pvt)),
                         (Float64(wat_mu_ref), Float64(co2_mu_ref)),
                         (lin(wat_mu_pvt), lin(co2_mu_pvt)),
                         dissolution, Float64(dis_rate), Float64(dis_max))
end

"""
    CO2LabVEFluid

Settings of co2lab's VE fluid; build with [`co2lab_ve_fluid`](@ref). Phase
order in all tuples: (brine, CO2). The type parameter is that of the
vertical model. `rho_pvt`/`mu_pvt` hold `(c, p_ref)` per
phase, `(0.0, 0.0)` for a constant property. `dissolution`, `dis_rate`,
`dis_max` as makeVEFluid's options.
"""
struct CO2LabVEFluid{M <: VerticalModel}
    vertical::M                               # vertical model (specification)
    res_water::Float64                        # residual brine saturation
    res_gas::Float64                          # residual CO2 saturation
    krw::Float64                              # brine relperm in the residual-CO2 zone
    krg::Float64                              # CO2 relperm in the plume
    rho_ref::NTuple{2, Float64}               # reference ("surface") densities
    rho_pvt::NTuple{2, NTuple{2, Float64}}    # (c, p_ref): ρ_ref·(1 + c·(p − p_ref))
    mu_ref::NTuple{2, Float64}                # reference viscosities
    mu_pvt::NTuple{2, NTuple{2, Float64}}     # (c, p_ref), as rho_pvt
    dissolution::Bool                         # CO2 dissolves into brine
    dis_rate::Float64                         # dissolution rate; 0 or Inf: instantaneous
    dis_max::Float64                          # saturated dissolved-CO2 ratio rs
end

"""
    instantaneous_dissolution(fluid) -> Bool

Whether `fluid` has dissolution with co2lab's instantaneous mode (`dis_rate`
0 or `Inf`).
"""
instantaneous_dissolution(f::CO2LabVEFluid) =
    f.dissolution && (f.dis_rate == 0 || isinf(f.dis_rate))

# `nothing` → constant; (c, p_ref) → linear; a 4-entry table spec is refused.
function _linear_pvt(spec)
    isnothing(spec) && return (0.0, 0.0)
    length(spec) == 2 ||
        error("co2lab's sampled PVT tables ($spec) are not supported; " *
              "give nothing (constant) or (c, p_ref) (linear)")
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
    _fine_saturations(m, i, zb, h, h_max, rw, rg) -> (s, smax)

CO2 saturation and maximum CO2 saturation of each fine cell of column `i`
(cumulative cell-bottom depths below the caprock `zb`) under vertical model
`m`. Returns two new vectors, top cell first.
"""
_fine_saturations(::SharpInterface, i, zb, h, h_max, rw, rg) =
    _column_saturations(zb, h, h_max, rw, rg)

"Linear property `ref·(1 + c·(p − p_ref))`; constant when `c = 0`."
_linear(ref, (c, p_ref), p) = ref * (1 + c * (p - p_ref))

"""
    _density_difference(fluid, p, rs) -> Real

Brine minus CO2 density at pressure `p` with dissolved-CO2 ratio `rs`, as
the model's black-oil densities (`update_bo_densities!`): the brine carries
its dissolved CO2.
"""
function _density_difference(fluid::CO2LabVEFluid, p, rs)
    (ρW, ρG), (pvt_w, pvt_g) = fluid.rho_ref, fluid.rho_pvt
    return _linear(1.0, pvt_w, p) * (ρW + rs * ρG) - _linear(1.0, pvt_g, p) * ρG
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

# The sharp interface does not read the densities, so it does not depend on
# them (an unchanged dependency graph and derivatives).
# Modifies: `hh` (plume heights).
@jutul_secondary function update_plume_heights!(hh, v::PlumeHeights{SharpInterface},
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

# The sharp interface: heights only, no saturations or densities.
# Modifies: `kr` (relative permeabilities).
@jutul_secondary function update_ve_relperm!(kr, v::VERelPerm{SharpInterface}, model,
                                            PlumeHeights, ColumnHeight, Rugosity, ix)
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
    _set_ve_fluid!(model, fluid)

Install co2lab's VE constitutive functions in the reservoir submodel `model`
(which must be on a domain from `ve_domain`).

Modifies: `model.parameters` (adds `:ColumnHeight`, `:Rugosity`),
`model.secondary_variables` (`:PlumeHeights`, `:RelativePermeabilities`,
`:CapillaryPressure`), `model.output_variables` (adds `:PlumeHeights`,
`:CapillaryPressure`). `fluid` is only read.
"""
function _set_ve_fluid!(model, fluid::CO2LabVEFluid)
    haskey(model.data_domain, :column_height) ||
        error("VE fluid needs domain[:column_height]; build the domain with ve_domain")

    prm, S = model.parameters, model.secondary_variables
    prm[:ColumnHeight] = ColumnHeight()
    prm[:Rugosity] = Rugosity()

    # per-cell VE flag; `nothing` (branch compiled away) when all cells are VE
    d = model.data_domain
    is_ve = haskey(d, :is_ve) && !all(d[:is_ve]) ? collect(Bool, d[:is_ve]) : nothing

    # the column-specific vertical model, shared by heights and relperm
    m = _vertical_model(fluid.vertical, d, fluid)
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

co2lab's PVT as Jutul variables: `b = 1 + c·(p − p_ref)` per phase (the
reference density is the "surface" density), and linear or constant
viscosities. Returns new variable objects; `fluid` is only read.
"""
function _co2lab_pvt(fluid::CO2LabVEFluid)
    (cw, pw), (cg, pg) = fluid.rho_pvt
    b = LinearShrinkage((1.0, 1.0), (cw, cg), (pw, pg))
    return (b, LinearViscosities(fluid.mu_ref, fluid.mu_pvt))
end
