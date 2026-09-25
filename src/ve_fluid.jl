# co2lab's VE constitutive functions and PVT (roadmap Stage 2 items 1 and 3).
#
# - plume heights (h, h_max) from sG, sGmax and the column height H
#   (co2lab's upscaledSat2height, homogeneous case);
# - the 'sharp_interface_simple' relative permeabilities, with endpoints and
#   the accretion layer of addVERelpermSharpInterface — from the caprock
#   rugosity, a property of the VE grid (`ve_domain`), per column;
# - pcWG = Δρ g h, added to the brine (reference-phase) pressure to give the
#   CO2 pressure; Δρ from our black-oil phase densities (researcher's
#   decision: natural in JutulDarcy; equals co2lab's without dissolved CO2);
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

# ---------------------------------- Settings ----------------------------------

"""
    co2lab_ve_fluid(; residual = (0.0, 0.0), krmax = (-1.0, -1.0),
                      co2_rho_ref = 760.0, wat_rho_ref = 1100.0,
                      co2_rho_pvt = nothing, wat_rho_pvt = nothing,
                      co2_mu_ref = 6e-5, wat_mu_ref = 8e-4,
                      co2_mu_pvt = nothing, wat_mu_pvt = nothing) -> CO2LabVEFluid

co2lab's VE fluid for the `'sharp_interface_simple'` relperm model, with
`makeVEFluid`'s option names and defaults, so an MRST case translates one to
one:
- `residual = (res_water, res_gas)`; `krmax = (krw, krg)` endpoint relperms,
  where a negative value means one minus the other phase's residual;
- `*_rho_ref`: reference densities, which are also the black-oil "surface"
  densities (so `b = ρ(p)/ρ_ref`);
- `*_rho_pvt`, `*_mu_pvt`: `nothing` (constant) or `(c, p_ref)` for
  `ref·(1 + c·(p − p_ref))`. co2lab's third option, a sampled EOS table (and
  its default for CO2 density), is not supported.

makeVEFluid's other options are not fluid properties and go elsewhere:
- `rugosity` → [`ve_domain`](@ref)`(...; rugosity)`, per column;
- `pvMult_p_ref`, `pvMult_fac` → [`setup_co2lab_model`](@ref)`(...;
  rock_compressibility = (reference_pressure = pvMult_p_ref,
  compressibility = pvMult_fac))`.
Passing them here throws an `ArgumentError` that says so.

Returns a new `CO2LabVEFluid`, with the endpoints resolved.
"""
function co2lab_ve_fluid(; residual = (0.0, 0.0), krmax = (-1.0, -1.0),
                           co2_rho_ref = 760.0, wat_rho_ref = 1100.0,
                           co2_rho_pvt = nothing, wat_rho_pvt = nothing,
                           co2_mu_ref = 6e-5, wat_mu_ref = 8e-4,
                           co2_mu_pvt = nothing, wat_mu_pvt = nothing,
                           rugosity = nothing, pvMult_p_ref = nothing,
                           pvMult_fac = nothing)
    _refuse_moved_options(rugosity, pvMult_p_ref, pvMult_fac)
    rw, rg = Float64.(residual)
    krw = krmax[1] < 0 ? 1 - rg : Float64(krmax[1])
    krg = krmax[2] < 0 ? 1 - rw : Float64(krmax[2])
    lin = _linear_pvt

    return CO2LabVEFluid(rw, rg, krw, krg,
                         (Float64(wat_rho_ref), Float64(co2_rho_ref)),
                         (lin(wat_rho_pvt), lin(co2_rho_pvt)),
                         (Float64(wat_mu_ref), Float64(co2_mu_ref)),
                         (lin(wat_mu_pvt), lin(co2_mu_pvt)))
end

# makeVEFluid options that belong to the grid or the model, not the fluid
function _refuse_moved_options(rugosity, pvMult_p_ref, pvMult_fac)
    isnothing(rugosity) ||
        throw(ArgumentError("rugosity is a property of the VE grid: give it to " *
                            "ve_domain(...; rugosity), per column or as a scalar"))
    isnothing(pvMult_p_ref) && isnothing(pvMult_fac) ||
        throw(ArgumentError("rock compressibility belongs to the model: give it to " *
                            "setup_co2lab_model(...; rock_compressibility = " *
                            "(reference_pressure = pvMult_p_ref, " *
                            "compressibility = pvMult_fac))"))
    return nothing
end

"""
    CO2LabVEFluid

Settings of co2lab's VE fluid; build with [`co2lab_ve_fluid`](@ref). Phase
order in all tuples: (brine, CO2). `rho_pvt`/`mu_pvt` hold `(c, p_ref)` per
phase, `(0.0, 0.0)` for a constant property.
"""
struct CO2LabVEFluid
    res_water::Float64
    res_gas::Float64
    krw::Float64
    krg::Float64
    rho_ref::NTuple{2, Float64}
    rho_pvt::NTuple{2, NTuple{2, Float64}}
    mu_ref::NTuple{2, Float64}
    mu_pvt::NTuple{2, NTuple{2, Float64}}
end

# `nothing` → constant; (c, p_ref) → linear; a 4-entry table spec is refused.
function _linear_pvt(spec)
    isnothing(spec) && return (0.0, 0.0)
    length(spec) == 2 ||
        error("co2lab's sampled PVT tables ($spec) are not supported; " *
              "give nothing (constant) or (c, p_ref) (linear)")
    return (Float64(spec[1]), Float64(spec[2]))
end

# ---------------------------------- Formulas ----------------------------------

"""
    _plume_heights(sG, sGmax, H, rw, rg) -> (h, h_max)

co2lab's `upscaledSat2height` for a vertically homogeneous column
(`upscaledSat2height.m`, `sharp_interface_case`): `h_max = sGmax·H/(1−rw)`,
`h = (H·sG − h_max·rg)/(1−rg−rw)`, clamped to `[0, H]`.
"""
function _plume_heights(sG, sGmax, H, rw, rg)
    h_max = sGmax * H / (1 - rw)
    h = (H * sG - h_max * rg) / (1 - rg - rw)
    return (min(max(h, zero(h)), H), h_max)
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
    krG = max(h - lost, zero(h)) * krg / H
    krW = (H - max(h_max, lost) + max(h_max - max(h, lost), zero(h)) * krw) / H
    return (krW, krG)
end

"Linear property `ref·(1 + c·(p − p_ref))`; constant when `c = 0`."
_linear(ref, (c, p_ref), p) = ref * (1 + c * (p - p_ref))

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

"Plume heights per column: row 1 `h`, row 2 `h_max` (m)."
struct PlumeHeights <: Jutul.VectorVariables
    res_water::Float64
    res_gas::Float64
end
Jutul.values_per_entity(model, ::PlumeHeights) = 2

# Uses sGmax_eff = max(min(sGmax, 1 − rw), sG): the stored sGmax is only
# updated after a converged step, so within a step the current sG must count
# too (MRST's updateState does this after every Newton update,
# CO2VEBlackOilTypeModel.m:352-353). Otherwise h_max < h while the plume grows.
# Modifies: `hh` (plume heights).
@jutul_secondary function update_plume_heights!(hh, v::PlumeHeights, model, sG, sGmax,
                                               ColumnHeight, ix)
    rw, rg = v.res_water, v.res_gas
    for i in ix
        smax = max(min(sGmax[i], 1 - rw), sG[i])
        hh[1, i], hh[2, i] = _plume_heights(sG[i], smax, ColumnHeight[i], rw, rg)
    end
end

"co2lab's `'sharp_interface_simple'` relperm; phases (brine, CO2)."
struct VESharpInterfaceRelPerm <: JutulDarcy.AbstractRelativePermeabilities
    krw::Float64
    krg::Float64
end

# Modifies: `kr` (relative permeabilities).
@jutul_secondary function update_ve_relperm!(kr, v::VESharpInterfaceRelPerm, model,
                                            PlumeHeights, ColumnHeight, Rugosity, ix)
    for i in ix
        h, h_max = PlumeHeights[1, i], PlumeHeights[2, i]
        kr[1, i], kr[2, i] = _sharp_interface_relperm(h, h_max, ColumnHeight[i],
                                                      v.krw, v.krg, Rugosity[i])
    end
end

"`pcWG = (ρ_brine − ρ_CO2)·g·h`, added to the brine pressure for CO2."
struct VECapillaryPressure <: JutulDarcy.AbstractCapillaryPressure end

# Modifies: `pc` (capillary pressure, one row: CO2 relative to brine).
@jutul_secondary function update_ve_pc!(pc, v::VECapillaryPressure, model, PlumeHeights,
                                       PhaseMassDensities, ix)
    g = Jutul.gravity_constant
    for i in ix
        pc[1, i] = (PhaseMassDensities[1, i] - PhaseMassDensities[2, i]) * g *
                   PlumeHeights[1, i]
    end
end

"Viscosities `μ_ref·(1 + c·(p − p_ref))` per phase (brine, CO2)."
struct LinearViscosities <: JutulDarcy.PhaseVariables
    mu_ref::NTuple{2, Float64}
    pvt::NTuple{2, NTuple{2, Float64}}
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

    S[:PlumeHeights] = PlumeHeights(fluid.res_water, fluid.res_gas)
    S[:RelativePermeabilities] = VESharpInterfaceRelPerm(fluid.krw, fluid.krg)
    S[:CapillaryPressure] = VECapillaryPressure()

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
