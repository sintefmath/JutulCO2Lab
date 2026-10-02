# Tabulated phase properties in pressure and temperature (roadmap Stage 7):
# JutulDarcy's CO2–brine tables (the default) and co2lab's sampled CO2 tables
# (for comparisons with co2lab).
#
# - JutulDarcy's: `co2_brine_property_tables` (its CO2Properties module), the
#   CSP11 tables by default (NIST-consistent; researcher, 2026-10-02), the Salo
#   et al. 2024 correlations on request and for salinity; density and
#   viscosity of both phases (pure brine, CO2);
# - co2lab's: the sampled tables of addSampledFluidProperties (CO2 by
#   makeVEFluid's default; water too where a case asks for it), read from an
#   MRST installation, interpolated bilinearly as co2lab's SampledProp2D does
#   (values clamped to the table). For comparisons with co2lab only: they
#   move to scripts/ when Stage 7 closes (researcher, 2026-10-02).
#
# A phase property is a `PhaseTable` (a function of p and T), or `(c, p_ref)`
# for co2lab's linear model `ref·(1 + c·(p − p_ref))`, `(0, 0)` for constant.

"""
    PhaseTable

A tabulated phase property `(p, T) -> value` (pressure in Pa, temperature in
K), from [`jutuldarcy_property_tables`](@ref) or
[`co2lab_sampled_table`](@ref).
"""
struct PhaseTable{F}
    f::F            # (p, T) -> property value
    source::Symbol  # where it comes from (:jutuldarcy, :co2lab)
end
(t::PhaseTable)(p, T) = t.f(p, T)

# --------------------------------- JutulDarcy ---------------------------------

"""
    jutuldarcy_property_tables(; co2_source = :csp11, co2_density = :nist,
                               salt_names = String[], salt_mole_fractions = Float64[])
        -> NamedTuple

JutulDarcy's CO2–brine property tables (`co2_brine_property_tables`), as
phase tables: `co2_source = :csp11` (default: the CSP11 tables, consistent
with NIST) or `:salo24` (the Salo et al. 2024 correlations, which take
salinity; their CO2 viscosity is 5–30 % above NIST's,
`analysis/stage7-table-properties-2026-10-02.txt`); `co2_density = :nist`
(default) or `:rk` (Redlich–Kwong). Salts (`salt_names`, e.g. `["NaCl"]`, and
their `salt_mole_fractions`) require `co2_source = :salo24`.

Returns the NamedTuple `(density, viscosity)`, each a tuple of new
[`PhaseTable`](@ref)s `(brine, CO2)`.
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

# ----------------------------------- co2lab -----------------------------------

"""
    co2lab_sampled_table(file) -> PhaseTable

One of co2lab's sampled property tables (a `.mat` file of
`addSampledFluidProperties`, e.g.
`CarbonDioxide_100000_400000000_278_524_800_800_D.mat` for CO2 density, `_V`
for viscosity; see [`co2lab_table_file`](@ref)), interpolated bilinearly on
its uniform (p, T) grid as co2lab's `SampledProp2D`, and clamped to its range
as makeVEFluid's tables are.

Returns a new [`PhaseTable`](@ref); the file is only read.
"""
function co2lab_sampled_table(file)
    d = MAT.matread(file)
    vals = Matrix{Float64}(d["vals"])
    pspan = Float64.(Tuple(d["P"]["span"]))
    tspan = Float64.(Tuple(d["T"]["span"]))
    return PhaseTable((p, T) -> _bilinear(vals, pspan, tspan, p, T), :co2lab)
end

"""
    co2lab_table_file(dir, property; fluid = :co2, p_range = (1e5, 4e8),
                      T_range = (278, 524), samples = (800, 800)) -> String

The path of one of co2lab's sampled tables in `dir` (e.g.
`mrst/core/examples/data/CO2lab_Sampled_Tables` of an MRST installation),
named as co2lab's `propFilename`: `property` `:density` or `:viscosity`,
`fluid` `:co2` or `:water`, the pressure range (Pa), temperature range (K)
and sample counts. The defaults are makeVEFluid's CO2 table. Returns the
path; it is not checked.
"""
function co2lab_table_file(dir, property::Symbol; fluid::Symbol = :co2,
                           p_range = (1e5, 4e8), T_range = (278, 524),
                           samples = (800, 800))
    tag = Dict(:density => "D", :viscosity => "V")[property]
    name = Dict(:co2 => "CarbonDioxide", :water => "Water")[fluid]
    # MATLAB's num2str(x, '%d_'): integers plainly, others in %e form
    num(x) = isinteger(x) ? string(Int(x)) : @sprintf("%.6e", x)
    parts = [name; num.(collect(p_range)); num.(collect(T_range));
             num.(collect(samples))]
    return joinpath(dir, join(parts, "_") * "_$tag.mat")
end

# co2lab's SampledProp2D: the cell of the uniform grid (the last one at the
# upper end), its local coordinates, the bilinear weights; clamped to the range
function _bilinear(vals, (p0, p1), (t0, t1), p, T)
    np, nt = size(vals) .- 1
    p = clamp(p, p0, p1)
    T = clamp(T, t0, t1)
    i = min(1 + floor(Int, Jutul.value((p - p0) * np / (p1 - p0))), np)
    j = min(1 + floor(Int, Jutul.value((T - t0) * nt / (t1 - t0))), nt)
    s = (p - (p0 + (i - 1) * (p1 - p0) / np)) / ((p1 - p0) / np)
    t = (T - (t0 + (j - 1) * (t1 - t0) / nt)) / ((t1 - t0) / nt)
    return (1 - s) * (1 - t) * vals[i, j] + s * (1 - t) * vals[i+1, j] +
           (1 - s) * t * vals[i, j+1] + s * t * vals[i+1, j+1]
end

# ------------------------------ Phase properties ------------------------------

"Phase property at pressure `p` and temperature `T`: tabulated, or linear."
_phase_value(ref, spec::PhaseTable, p, T) = spec(p, T)
_phase_value(ref, (c, p_ref)::NTuple{2, Float64}, p, T) = ref * (1 + c * (p - p_ref))

"Whether a fluid needs a temperature (any tabulated property)."
_tabulated(specs...) = any(s -> s isa PhaseTable, specs)

# ------------------------------- Jutul variables -------------------------------

"""
Shrinkage factors `b = ρ(p, T)/ρ_ref` per phase (brine, CO2), each tabulated
or linear (`pvt`), at the cell's pressure and temperature (the `Temperature`
parameter).
"""
struct PVTShrinkage{R} <: JutulDarcy.PhaseVariables
    rho_ref::NTuple{2, Float64}  # reference ("surface") densities
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
Viscosities per phase (brine, CO2), each tabulated or linear (`pvt`, with
`mu_ref` for linear ones), at the cell's temperature and, with `P = true`, at
the phase pressures as co2lab's (CO2 at `p + pcWG`; wells, `P = false`, at
`p`).
"""
struct PVTViscosities{U, P} <: JutulDarcy.PhaseVariables
    mu_ref::NTuple{2, Float64}  # reference viscosities (linear ones)
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
