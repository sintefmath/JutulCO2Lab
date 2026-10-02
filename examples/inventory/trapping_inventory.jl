# Where does the injected CO2 end up? The trapping inventory of a VE simulation
# (co2lab-ve's postprocessStates and plotTrappingDistribution): CO2 injected
# into the Johansen formation of the CO2 Storage Atlas for 20 years, followed
# for 180 more, split by trapping mechanism — structural traps (from the
# spill-point analysis), residual trapping (hysteresis), dissolution.
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/inventory/trapping_inventory.jl
#
# Plotting backend: GLMakie (a window per figure, 3D views rotatable with the
# mouse; close it to continue) by default; JUTULCO2LAB_BACKEND=CairoMakie
# writes the figures to files instead.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Jutul, JutulDarcy, Printf

# ----------------------- 1. The formation and its traps -----------------------

# Johansen, coarsened by 4 (2 km columns), as a VE domain with the Atlas's
# rock. A rugosity (atlas_ve_domain(...; rugosity)) would add subscale traps;
# here the caprock is smooth below the grid's resolution.
case = atlas_ve_domain("Johansenfm"; coarsening = 4)
ve = case.ve

# Structural traps of the caprock, node-based as in co2lab's VE examples.
s = trap_surface(case.mesh, case.columns)
traps = trap_analysis(s; method = :node)
volumes = trap_volumes(s, traps; poro = ve.domain[:porosity])
@printf("%d columns, %d traps holding %.2e m3 of pore space\n", length(s.zc),
        length(volumes), sum(volumes))

# ----------------------------- 2. Fluid and model -----------------------------

# Brine and CO2 with residual saturations (hysteresis: CO2 left behind by
# the retreating plume is trapped), linear densities, and dissolution of CO2
# into the brine at a finite rate.
g, rho_w, rho_co2 = 9.80665, 1020.0, 700.0
fluid = co2lab_ve_fluid(default_pvt = nothing,
                        residual = (0.27, 0.20), krmax = (1.0, 1.0),
                        co2_rho_ref = rho_co2, wat_rho_ref = rho_w,
                        co2_rho_pvt = (4.5e-9, 100e5), wat_rho_pvt = (4.6e-10, 100e5),
                        co2_mu_ref = 6e-5, wat_mu_ref = 8e-4,
                        dissolution = true, dis_rate = 5e-11, dis_max = 0.03)

# One well at the Johansen benchmark's injection point, 1 Mt of CO2 a year.
w = johansen_well_point()
well_column = argmin((s.cx .- w[1]) .^ 2 .+ (s.cy .- w[2]) .^ 2)
well = setup_ve_well(ve, well_column; name = :inj, radius = 0.1)
rate = 1e9 / (365.25 * 86400) / rho_co2
inject = co2_injector_control(TotalRateTarget(rate), fluid)

# All outer sides open at hydrostatic pressure.
sides, _ = ve_boundary_sides(case.mesh, case.columns, ve.domain)
bc = ve_boundary_conditions(ve, sides, ve_select_sides(sides, case.columns);
                            pressure = z -> rho_w * g * z)
model = setup_co2lab_model(ve.domain; wells = [well], fluid)

# -------------------------------- 3. Simulate ---------------------------------

# 20 years of injection in yearly steps, then 180 years in 10-year steps.
year = 365.25 * 86400
f_inj = setup_reservoir_forces(model; control = Dict(:inj => inject), bc = bc)
f_shut = setup_reservoir_forces(model; control = Dict(:inj => DisabledControl()),
                                bc = bc)
dt = [fill(year, 20); fill(10year, 18)]
forces = [k <= 20 ? f_inj : f_shut for k in eachindex(dt)]
p0 = ve_hydrostatic_pressure(ve.geometry; rho_w = rho_w, g = g)
state0 = setup_co2lab_state(model; pressure = p0, forces = f_inj)
prm = setup_parameters(model)
run = JutulCase(model, dt, forces; state0 = state0, parameters = prm)
sim, cfg = setup_co2lab_simulator(run; info_level = -1)
result = simulate_reservoir(run; simulator = sim, config = cfg)

# ------------------------- 4. The trapping inventory --------------------------

# Per report step: dissolved, residual (in and out of traps, and the part of
# the plume that will stay behind as residual), structural (plume in traps),
# free plume, and what has left the formation.
history = trapping_history(result, model, fluid, s, traps)
for k in (5, 20, 30, length(history.t))
    @printf("\nyear %.0f (Mt):\n", history.t[k] / year)
    for c in TRAPPING_CATEGORIES
        @printf("  %-20s %7.3f\n", c, history.inventory[k][c] / 1e9)
    end
    @printf("  %-20s %7.3f\n", :exited, history.exited[k] / 1e9)
end
fig = plot_trapping_history(history; title = "Johansen: trapping inventory")
show_figure(fig, "johansen_inventory")

# ----------------------- 5. The plume at the end, in 3D -----------------------

# The traps and their spill regions on the caprock, and over them the columns
# with mobile CO2 at the end, coloured by plume height; the well as a line.
h = result.states[end][:PlumeHeights][1, :]
plume = findall(>(0.01), h)
fig = Figure(size = (1200, 900))
years = round(Int, history.t[end] / year)
ax = ve_axis3(fig[1, 1], s; view = (-62, 60), exaggeration = 5,
              title = "traps, and the plume after $years years")
plot_traps!(ax, s, traps; rivers = false)
# the plume a few metres above the caprock, so that it is drawn over it
plot_ve_surface!(ax, (; s..., polygons = s.polygons[plume]); color = h[plume],
                 colormap = :plasma, offset = -3.0, shading = true, strokewidth = 0)
plot_ve_map!(ax, s; wells = [well_column], outline = false)
Colorbar(fig[1, 2]; limits = extrema(h[plume]), colormap = :plasma,
         label = "plume height (m)")
show_figure(fig, "johansen_final_plume")
