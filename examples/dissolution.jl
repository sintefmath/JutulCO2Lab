# co2lab's sloping aquifer with residual trapping and CO2 dissolving into the
# brine (sloping_aquifer.m's 'dissolution' options), simulated with JutulCO2Lab
# three times: without dissolution, with co2lab's rate-driven dissolution, and
# with instantaneous dissolution. The examples sloping_aquifer.jl and
# residual_trapping.jl explain the setup.
#
#   julia --project=. examples/dissolution.jl
#
# No input files: co2lab's case is built by `sloping_aquifer` (40 layers); CO2
# properties from JutulDarcy's tables at its reference conditions.

using JutulCO2Lab, Jutul, JutulDarcy, Printf

# ------------------------------- 1. The case ----------------------------------

aq = sloping_aquifer(; nz = 40)
g = Jutul.gravity_constant
rho_w, rho_co2 = aq.rho_w, aq.co2_rho           # brine; CO2 at reference p and T

mesh, columns, ve = aq.mesh, aq.columns, aq.ve
geo = ve.geometry

# ------------------------------- 2. Three fluids -------------------------------

# Residual saturations as in residual_trapping.jl (hysteresis on).
srw, srg = 0.25, 0.25

# Brine holds at most dis_max m³ of CO2 per m³ of brine, both at reference
# conditions (black-oil's rs).
dis_max = 0.05

# Rate-driven: a column's brine takes up CO2 at dis_rate · pore volume / column
# height while the column has free CO2 (co2lab's rate constant).
# Instantaneous (dis_rate = 0): a column's brine is saturated as soon as it
# meets CO2, as far as the CO2 there suffices.
modes = (without = nothing, rate_driven = 3e-11, instantaneous = 0.0)

function make_fluid(dis_rate)
    dis = isnothing(dis_rate) ? (;) :
          (dissolution = true, dis_rate = dis_rate, dis_max = dis_max)
    return co2lab_ve_fluid(; default_pvt = nothing,
                           residual = (srw, srg), krmax = (1.0, (1 - srw)^2),
                           co2_rho_ref = rho_co2, wat_rho_ref = rho_w,
                           co2_mu_ref = aq.co2_mu, wat_mu_ref = aq.mu_w, dis...)
end

rock_c = (reference_pressure = aq.p_ref, compressibility = aq.rock_compressibility)

# -------------------------- 3. Well, boundary, model ---------------------------

well_column = aq.well_column
well = setup_ve_well(ve, well_column; name = :inj, radius = aq.well_radius)

# The updip end is open, at hydrostatic brine pressure. Brine leaving there
# takes its dissolved CO2 along (the default, `dissolved = :carry`).
sides, open_sides = aq.sides, aq.open_sides

# ------------------------------- 4. Schedule ----------------------------------

# co2lab's 2 years of injection and 2 of migration, then 46 more years of
# migration in 2-year steps: dissolution is slow.
year = si_unit(:year)
dt = vcat(aq.dt, fill(2year, 23))
injecting = vcat(aq.injecting, fill(false, 23))
years = cumsum(dt) ./ year

# ---------------------------------- 5. Run ------------------------------------

# One run: model, forces, initial state (wellbore filled with the injected CO2),
# simulation. Returns the reservoir states and the Newton iteration count.
function run(fluid)
    bc = ve_boundary_conditions(ve, sides, open_sides; pressure = z -> rho_w * g * z)
    model = setup_co2lab_model(ve.domain; wells = [well], fluid = fluid,
                               rock_compressibility = rock_c)

    inject = co2_injector_control(TotalRateTarget(aq.mass_rate / rho_co2), fluid)
    on = setup_reservoir_forces(model; control = Dict(:inj => inject), bc = bc)
    off = setup_reservoir_forces(model; control = Dict(:inj => DisabledControl()),
                                 bc = bc)
    forces = [i ? on : off for i in injecting]

    p0 = ve_hydrostatic_pressure(geo; rho_w = rho_w, g = g)
    state0 = setup_co2lab_state(model; pressure = p0, forces = on)

    case = JutulCase(model, dt, forces; state0 = state0,
                     parameters = setup_parameters(model))
    sim, cfg = setup_co2lab_simulator(case; info_level = -1)
    result = simulate_reservoir(case; simulator = sim, config = cfg).result

    newton = count(haskey(s, :linear_iterations) for r in result.reports
                   for m in r[:ministeps] for s in m[:steps])
    return (states = [s[:Reservoir] for s in result.states], newton = newton)
end

runs = map(dis_rate -> (@time run(make_fluid(dis_rate))), modes)
println("Newton iterations: ", join(("$k $(r.newton)" for (k, r) in pairs(runs)), ", "))

# ----------------------- 6. Where the CO2 is, over time -----------------------

# Per column: the free CO2 splits into the mobile plume (height h) and the
# residual zone below it (h to h_max), as in residual_trapping.jl; the
# dissolved CO2 is the state's DissolvedCO2 (no dissolution: none).
function inventory(s)
    h, h_max = s[:PlumeHeights][1, :], s[:PlumeHeights][2, :]
    dissolved = haskey(s, :DissolvedCO2) ? s[:DissolvedCO2] : zero(h)
    # TotalMasses holds all CO2 in a column, free and dissolved
    free = s[:TotalMasses][2, :] .- dissolved
    mobile_vol, trapped_vol = (1 - srw) .* h, srg .* max.(h_max .- h, 0.0)
    share = trapped_vol ./ max.(mobile_vol .+ trapped_vol, eps())
    return (mobile = sum(free .* (1 .- share)), residual = sum(free .* share),
            dissolved = sum(dissolved))
end

# exited: injected minus what is in place (−0.00: rounding of a tiny balance)
injected = cumsum([i ? aq.mass_rate * t : 0.0
                   for (i, t) in zip(injecting, dt)])
report_steps = [argmin(abs.(years .- y)) for y in (2, 4, 10, 20, 50)]

for (mode, r) in pairs(runs)
    println("\n$(replace(string(mode), '_' => '-')) dissolution")
    println("  year   injected     mobile   residual  dissolved     exited   (kt)")
    for k in report_steps
        c = inventory(r.states[k])
        exited = injected[k] - c.mobile - c.residual - c.dissolved
        @printf("  %4.0f   %8.2f   %8.2f   %8.2f   %8.2f   %8.2f\n", years[k],
                injected[k] / 1e6, c.mobile / 1e6, c.residual / 1e6,
                c.dissolved / 1e6, exited / 1e6)
    end
end

# ----------------------- 7. The plume after 50 years --------------------------

# Dissolution eats into the residual zone first; instantaneous dissolution also
# shortens the plume. Brine carrying dissolved CO2 spreads beyond the plume.
println()
for (mode, r) in pairs(runs)
    h = r.states[end][:PlumeHeights][1, :]
    rs = get(r.states[end], :Rs, zero(h))
    @printf("%-14s mobile CO2 in %2d columns; ", mode, count(>(1e-3), h))
    @printf("dissolved CO2 in %2d columns, saturated in %2d\n",
            count(>(1e-6), rs), count(>=(dis_max - 1e-6), rs))
end

# --------------------------- 8. Back to the 3D grid ---------------------------

# The 3D picture of the rate-driven run: dissolved CO2 per fine cell.
fluid = make_fluid(modes.rate_driven)
fine = ve_reconstruct(runs.rate_driven.states[end], ve, columns, fluid)
@printf("rate-driven, 3D cells: %d hold free CO2, %d hold dissolved CO2 (of %d)\n",
        count(>(1e-6), fine.sG), count(>(1e-6), fine.rs), length(fine.sG))
