# Phase properties (density, viscosity) as functions of pressure and
# temperature, for the VE and fine-scale fluids.
#
# A phase property is either
# - a `PhaseTable`, a function of `(p, T)`: JutulDarcy's CO2–brine tables
#   (`jutuldarcy_property_tables`, the default), or co2lab's sampled tables
#   (read by the comparison scripts, for comparisons with co2lab only); or
# - `(c, p_ref)`, co2lab's linear model `ref·(1 + c·(p − p_ref))`; `(0, 0)`
#   for a constant.

"""
    PhaseTable

A tabulated phase property, callable as `t(p, T)` (pressure in Pa,
temperature in K). Made by [`jutuldarcy_property_tables`](@ref), or wraps
any function of `p` and `T` (e.g. co2lab's sampled tables).

# Fields

- `f`: the property as a function `(p, T) -> value`.
- `source`: a label for where the table comes from (e.g. `:jutuldarcy`,
  `:co2lab`).
"""
struct PhaseTable{F}
    f::F            # (p, T) -> property value
    source::Symbol  # where it comes from (e.g. :jutuldarcy)
end
(t::PhaseTable)(p, T) = t.f(p, T)

# --------------------------------- JutulDarcy ---------------------------------

"""
    jutuldarcy_property_tables(; co2_source = :csp11, co2_density = :nist,
                               salt_names = String[], salt_mole_fractions = Float64[])
        -> NamedTuple

JutulDarcy's CO2–brine property tables (`co2_brine_property_tables`) as
[`PhaseTable`](@ref)s, for tabulated densities and viscosities of pure
brine and CO2.

# Arguments

- `co2_source`: `:csp11` (default; the CSP11 tables, consistent with NIST)
  or `:salo24` (the Salo et al. 2024 correlations; their CO2 viscosity is
  5–30 % above NIST's).
- `co2_density`: `:nist` (default) or `:rk` (Redlich–Kwong).
- `salt_names`, `salt_mole_fractions`: salts (e.g. `["NaCl"]`) and their
  mole fractions; need `co2_source = :salo24` (an error otherwise).

# Returns

The NamedTuple `(density, viscosity)`, each a tuple of new `PhaseTable`s
`(brine, CO2)`: density in kg/m³, viscosity in Pa·s.
"""
function jutuldarcy_property_tables(; co2_source = :csp11, co2_density = :nist,
                                    salt_names = String[],
                                    salt_mole_fractions = Float64[])
    # JutulDarcy ignores salts with a table (only a message): refuse instead
    isempty(salt_names) || co2_source == :salo24 ||
        error("salinity needs co2_source = :salo24 (the CSP11 tables are for " *
              "pure water)")
    tab = JutulDarcy.CO2Properties.co2_brine_property_tables(; co2_source,
                                                             co2_density, salt_names,
                                                             salt_mole_fractions)
    phase(t, ph) = PhaseTable((p, T) -> t(p, T)[ph], :jutuldarcy)
    rho, mu = tab[:density], tab[:viscosity]
    return (density = (phase(rho, 1), phase(rho, 2)),
            viscosity = (phase(mu, 1), phase(mu, 2)))
end

# ------------------------------ Phase properties ------------------------------

"Phase property at pressure `p` and temperature `T`: tabulated, or linear."
_phase_value(ref, spec::PhaseTable, p, T) = spec(p, T)
_phase_value(ref, (c, p_ref)::NTuple{2, Float64}, p, T) = ref * (1 + c * (p - p_ref))

"Whether a fluid needs a temperature (any tabulated property)."
_tabulated(specs...) = any(s -> s isa PhaseTable, specs)

# ------------------------------- Jutul variables -------------------------------

"""
    PVTShrinkage

Shrinkage factors `b = ρ(p, T)/ρ_ref` per phase (brine, CO2), at the cell's
pressure and `Temperature` parameter; each phase tabulated or linear.

# Fields

- `rho_ref`: reference ("surface") densities `ρ_ref` (brine, CO2) (kg/m³);
  also the reference value of a linear density model.
- `pvt`: the density model per phase (brine, CO2): a [`PhaseTable`](@ref),
  or `(c, p_ref)` for co2lab's linear model `ρ_ref·(1 + c·(p − p_ref))`.
"""
struct PVTShrinkage{R} <: JutulDarcy.PhaseVariables
    rho_ref::NTuple{2, Float64}  # reference ("surface") densities (kg/m³)
    pvt::R                       # per phase: PhaseTable or (c, p_ref)
end

# Modifies: `b` (shrinkage factors).
@jutul_secondary function update_pvt_shrinkage!(b, v::PVTShrinkage, model, Pressure,
                                               Temperature, ix)
    for i in ix, ph in 1:2
        b[ph, i] = _phase_value(v.rho_ref[ph], v.pvt[ph], Pressure[i],
                                Temperature[i]) / v.rho_ref[ph]
    end
end

"""
    PVTViscosities(mu_ref, pvt; phase_pressures = true)

Viscosities per phase (brine, CO2), at the cell's `Temperature` parameter;
each phase tabulated or linear (`mu_ref` is used for linear ones only).

The pressure: with `phase_pressures = true` (type parameter `P`; the
reservoir) each phase's own, as co2lab (CO2 at `p + pcWG`); with `false`
(wells, which have no capillary pressure) `p`.

# Fields

- `mu_ref`: reference viscosities (brine, CO2) (Pa·s), used by linear phases
  only.
- `pvt`: the viscosity model per phase (brine, CO2): a [`PhaseTable`](@ref),
  or `(c, p_ref)` for co2lab's linear model `mu_ref·(1 + c·(p − p_ref))`.
"""
struct PVTViscosities{U, P} <: JutulDarcy.PhaseVariables
    mu_ref::NTuple{2, Float64}  # reference viscosities (Pa·s; linear ones)
    pvt::U                      # per phase: PhaseTable or (c, p_ref)
end
PVTViscosities(mu_ref, pvt::U; phase_pressures::Bool = true) where {U} =
    PVTViscosities{U, phase_pressures}(mu_ref, pvt)

# Modifies: `mu` (phase viscosities).
@jutul_secondary function update_pvt_viscosities!(mu, v::PVTViscosities{<:Any, true},
                                                 model, Pressure, CapillaryPressure,
                                                 Temperature, ix)
    for i in ix
        p, T = Pressure[i], Temperature[i]
        mu[1, i] = _phase_value(v.mu_ref[1], v.pvt[1], p, T)
        mu[2, i] = _phase_value(v.mu_ref[2], v.pvt[2], p + CapillaryPressure[1, i], T)
    end
end

# Modifies: `mu` (phase viscosities).
@jutul_secondary function update_pvt_viscosities!(mu, v::PVTViscosities{<:Any, false},
                                                 model, Pressure, Temperature, ix)
    for i in ix, ph in 1:2
        mu[ph, i] = _phase_value(v.mu_ref[ph], v.pvt[ph], Pressure[i], Temperature[i])
    end
end
