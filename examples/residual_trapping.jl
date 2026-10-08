# co2lab's sloping aquifer with residual trapping (slopingAquiferResidualCapillary.m,
# its sharp-interface configuration: residual saturations 0.25/0.25, zres = 40),
# simulated with JutulCO2Lab: one run, top to bottom. The Stage 2 example
# (sloping_aquifer.jl) explains the setup; this one adds hysteresis.
#
#   scripts/julia --project=JutulCO2Lab JutulCO2Lab/examples/residual_trapping.jl
#
# No input files: co2lab's case is built by `sloping_aquifer` (40 layers); CO2
# properties from JutulDarcy's tables at its reference conditions.

using JutulCO2Lab, Jutul, JutulDarcy, Printf

# ------------------------------- 1. The case ----------------------------------

aq = sloping_aquifer(; nz = 40)
g = Jutul.gravity_constant
rho_w, rho_co2 = aq.rho_w, aq.co2_rho           # brine; CO2 at reference p and T

# --------------------------------- 2. VE grid ---------------------------------

mesh, columns, ve = aq.mesh, aq.columns, aq.ve
geo = ve.geometry

# --------------------- 3. Fluid with residual saturations ---------------------

# Brine stays behind where CO2 has been (srw); CO2 stays behind where brine
# returns (srg). srg > 0 switches hysteresis on: sGmax becomes an unknown.
srw, srg = 0.25, 0.25

# CO2 relperm endpoint: the example's quadratic curve at 1 - srw.
fluid = co2lab_ve_fluid(default_pvt = nothing,
                        residual = (srw, srg), krmax = (1.0, (1 - srw)^2),
                        co2_rho_ref = rho_co2, wat_rho_ref = rho_w,
                        co2_mu_ref = aq.co2_mu, wat_mu_ref = aq.mu_w)

rock_c = (reference_pressure = aq.p_ref, compressibility = aq.rock_compressibility)

# ----------------------- 4. Well, open boundary, model ------------------------

well_column = aq.well_column
well = setup_ve_well(ve, well_column; name = :inj, radius = aq.well_radius)
inject = co2_injector_control(TotalRateTarget(aq.mass_rate / rho_co2), fluid)

# The updip end (x = L) is open, at hydrostatic brine pressure.
sides, open_sides = aq.sides, aq.open_sides
bc = ve_boundary_conditions(ve, sides, open_sides; pressure = z -> rho_w * g * z)

model = setup_co2lab_model(ve.domain; wells = [well], fluid = fluid,
                           rock_compressibility = rock_c)

# The hysteresis unknown and its equation.
r = model[:Reservoir]
println("unknowns: ", join(keys(r.primary_variables), ", "),
        ";  extra equation: ", r.equations[:hysteresis] isa HysteresisEquation)

# ----------------------- 5. Schedule and initial state ------------------------

f_inject = setup_reservoir_forces(model; control = Dict(:inj => inject), bc = bc)
f_shut = setup_reservoir_forces(model; control = Dict(:inj => DisabledControl()),
                                bc = bc)

dt = aq.dt                                        # 31 injection + 12 migration steps
injecting = aq.injecting
forces = [i ? f_inject : f_shut for i in injecting]

# Brine at hydrostatic pressure, no CO2; sGmax starts at sG = 0.
p0 = ve_hydrostatic_pressure(geo; rho_w = rho_w, g = g)
state0 = setup_co2lab_state(model; pressure = p0, forces = f_inject)

# -------------------------------- 6. Simulate ---------------------------------

prm = setup_parameters(model)
case = JutulCase(model, dt, forces; state0 = state0, parameters = prm)
sim, cfg = setup_co2lab_simulator(case; info_level = -1)
result = simulate_reservoir(case; simulator = sim, config = cfg).result
states = [s[:Reservoir] for s in result.states]

# ----------------------- 7. Where the CO2 is, over time -----------------------

# Per column, the plume (height h, CO2 at 1 - srw) is mobile; the residual zone
# below it (from h down to h_max, CO2 at srg) is trapped. Split each column's
# CO2 mass in that ratio.
function split_co2(s)
    h, h_max = s[:PlumeHeights][1, :], s[:PlumeHeights][2, :]
    mobile_vol = (1 - srw) .* h
    trapped_vol = srg .* max.(h_max .- h, 0.0)
    m = s[:TotalMasses][2, :]
    share = trapped_vol ./ max.(mobile_vol .+ trapped_vol, eps())
    return sum(m .* (1 .- share)), sum(m .* share)
end

injected = cumsum([i ? aq.mass_rate * t : 0.0
                   for (i, t) in zip(injecting, dt)])
years = cumsum(dt) ./ (365.25 * 86400)

println("\n  year   injected     mobile    trapped     exited   (kt)")
for k in (10, 31, 34, 37, 40, 43)
    mobile, trapped = split_co2(states[k])
    exited = injected[k] - mobile - trapped
    @printf("  %4.1f   %8.2f   %8.2f   %8.2f   %8.2f\n", years[k],
            injected[k] / 1e6, mobile / 1e6, trapped / 1e6, exited / 1e6)
end

# -------------------------- 8. The trail left behind --------------------------

final = states[end]
h, h_max = final[:PlumeHeights][1, :], final[:PlumeHeights][2, :]

# Columns the plume has passed and left: residual CO2 only.
trail = findall(i -> h[i] < 1e-3 && h_max[i] > 1e-3, eachindex(h))
plume = findall(>(1e-3), h)
a, b = extrema(plume)
@printf("\nplume now in columns %d-%d (mobile CO2 in %d of them)\n", a, b, length(plume))
@printf("residual CO2 only: %d columns downslope, %d between plume parts, %d updip\n",
        count(<(a), trail), count(i -> a < i < b, trail), count(>(b), trail))

# In the trail, sG = srg * h_max / H, and sGmax remembers the plume that was there.
i = trail[1]
@printf("column %d: sG %.4f = srg * h_max / H = %.4f;  sGmax %.4f\n",
        i, final[:sG][i], srg * h_max[i] / geo.H[i], final[:sGmax][i])

# --------------------------- 9. Back to the 3D grid ---------------------------

# The 3D picture: the plume under the caprock (CO2 at 1 - srw), and below it
# the residual zone (CO2 at srg); sGmax per 3D cell shows where CO2 has been.
fine = ve_reconstruct(final, ve, columns, fluid)
holding = count(>(1e-6), fine.sG)
in_plume = count(>(srg + 1e-6), fine.sG)
@printf("3D cells holding CO2: %d (%d reach into the plume, %d residual only)\n",
        holding, in_plume, holding - in_plume)
@printf("3D cells CO2 has reached (sGmax > 0): %d of %d\n",
        count(>(1e-6), fine.sGmax), length(fine.sGmax))
