# co2lab's sloping aquifer (slopingAquiferSharpInterface.m, cross-section, zres = 5),
# simulated with JutulCO2Lab: one run, top to bottom.
#
#   scripts/julia --project=JutulCO2Lab JutulCO2Lab/examples/sloping_aquifer.jl
#
# No input files: the case is built here as co2lab's sloping_aquifer.m builds
# it (its grid, make_testgrid, matches co2lab's node for node).

using JutulCO2Lab, Jutul, JutulDarcy, Printf

# ------------------------------- 1. The case ----------------------------------

# co2lab's case: grid, rock, well position, open boundary, schedule and
# constants. CO2's density and viscosity at the case's reference pressure and
# 70 °C come from JutulDarcy's tables (co2lab's own differ by < 0.2 %).
aq = sloping_aquifer(; nz = 5)
g = Jutul.gravity_constant
rho_w, rho_co2 = aq.rho_w, aq.co2_rho           # brine; CO2 at reference p and T

# --------------------------------- 2. VE grid ---------------------------------

# The 3D grid, a Jutul mesh built from a corner-point grid section.
mesh = aq.mesh

# Its 3D cells grouped into vertical columns; each column becomes one VE cell.
columns = aq.columns

# The VE domain: one coarse cell per column, thickness-averaged rock, gravity
# taken along the caprock. `geo` has caprock depth z_top and height H per column.
ve = aq.ve
geo = ve.geometry

# ---------- 3. Fluid (sharp interface, constant densities) and rock -----------

fluid = co2lab_ve_fluid(default_pvt = nothing, krmax = (1.0, 1.0),
                        co2_rho_ref = rho_co2, wat_rho_ref = rho_w,
                        co2_mu_ref = aq.co2_mu, wat_mu_ref = aq.mu_w)

# Rock compressibility belongs to the model (co2lab puts it in the fluid).
rock_c = (reference_pressure = aq.p_ref, compressibility = aq.rock_compressibility)

# ----------------------------- 4. Injection well ------------------------------

# One VE well in co2lab's well column, open over the whole column height.
well_column = aq.well_column
well = setup_ve_well(ve, well_column; name = :inj, radius = aq.well_radius)

# Rate given as volume at the CO2 reference density (co2lab's "surface" volume).
inject = co2_injector_control(TotalRateTarget(aq.mass_rate / rho_co2), fluid)

# ------------------------------ 5. Open boundary ------------------------------

# co2lab opens only the updip end of the cross-section (x = L: the :right side,
# as `ve_select_sides(sides, columns; directions = (:right,))` picks it).
sides, open_sides = aq.sides, aq.open_sides

# Hydrostatic brine pressure at the side's caprock edge; brine flows in if needed.
bc = ve_boundary_conditions(ve, sides, open_sides; pressure = z -> rho_w * g * z)

# -------------------- 6. Model, schedule and initial state --------------------

model = setup_co2lab_model(ve.domain; wells = [well], fluid = fluid,
                           rock_compressibility = rock_c)

# Also store what the boundary fluxes are computed from (for step 8).
push!(model[:Reservoir].output_variables, :RelativePermeabilities,
      :PhaseViscosities, :PhaseMassDensities)

# Two controls: injecting, then shut in; the boundary stays open throughout.
f_inject = setup_reservoir_forces(model; control = Dict(:inj => inject), bc = bc)
shut = DisabledControl()
f_shut = setup_reservoir_forces(model; control = Dict(:inj => shut), bc = bc)

dt = aq.dt                                        # 31 injection + 12 migration steps
injecting = aq.injecting
forces = [i ? f_inject : f_shut for i in injecting]

# Brine at hydrostatic pressure, no CO2 (pressure is at the caprock); the
# injector's wellbore starts filled with the CO2 it injects.
p0 = ve_hydrostatic_pressure(geo; rho_w = rho_w, g = g)
state0 = setup_co2lab_state(model; pressure = p0, forces = f_inject)

# -------------------------------- 7. Simulate ---------------------------------

prm = setup_parameters(model)
case = JutulCase(model, dt, forces; state0 = state0, parameters = prm)

# Keep the simulator's sub-steps: needed to integrate the boundary outflow exactly.
sim, cfg = setup_co2lab_simulator(case; info_level = 0, output_substates = true)
result = simulate_reservoir(case; simulator = sim, config = cfg).result

# ------------------------- 8. Results on the VE grid --------------------------

# CO2 injected and exited through the boundary: rates × sub-step lengths.
substates, dts, _ = Jutul.expand_to_ministeps(result)
injected = sum(s[:Facility][:TotalSurfaceMassRate][1] * t
               for (s, t) in zip(substates, dts))
exited = sum(ve_boundary_mass_fluxes(bc, s[:Reservoir])[2, 1] * t
             for (s, t) in zip(substates, dts))

# CO2 still in the aquifer at the end. The wellbore starts and ends filled
# with CO2, so it is not part of the budget.
final = result.states[end]
in_place = sum(final[:Reservoir][:TotalMasses][2, :])
in_well = final[:inj][:TotalMasses][2, 1]

@printf("CO2 injected           %10.4e kg\n", injected)
@printf("  in the aquifer       %10.4e kg\n", in_place)
@printf("  exited (updip end)   %10.4e kg\n", exited)
@printf("  unaccounted          %10.1e (relative)\n", (in_place + exited) / injected - 1)
@printf("(the wellbore holds %.0f kg of CO2 throughout)\n", in_well)

# The plume at the end: saturation and CO2 height under the caprock, per column.
res = final[:Reservoir]
plume = findall(>(1e-3), res[:sG])
@printf("plume in columns %d-%d; thickest %.2f m (column height %.0f m)\n",
        extrema(plume)..., maximum(res[:PlumeHeights][1, :]), geo.H[1])

# Pressure: back to hydrostatic after injection stops?
@printf("largest pressure change from initial: %.4f bar\n",
        maximum(abs.(res[:Pressure] .- p0)) / 1e5)

# The well's bottom-hole pressure when injection stops.
last_injecting = findlast(injecting)
bhp = result.states[last_injecting][:Facility][:BottomHolePressure][1]
@printf("injector BHP at end of injection: %.2f bar\n", bhp / 1e5)

# --------------- 9. Back to the 3D grid (co2lab's VEstates23D) ----------------

# Saturation per 3D cell from the plume height; pressure hydrostatic from the caprock.
fine = ve_reconstruct(res, ve, columns, fluid)
@printf("3D cells holding CO2: %d of %d (the top layers under the plume)\n",
        count(>(1e-3), fine.sG), length(fine.sG))
@printf("3D pressure range: %.2f - %.2f bar\n", extrema(fine.pressure) ./ 1e5...)
