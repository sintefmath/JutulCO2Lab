# co2lab's exampleVEJohansen.m with JutulCO2Lab: CO2 injected into the
# Johansen formation (the NPD sector model) for 50 years and followed for 950
# more, in a vertical-equilibrium model with a capillary fringe — the CO2 and
# brine meet in a transition zone set by the rock's capillary pressure, not at
# a sharp interface — and residual trapping. At the end, where the CO2 is, by
# trapping mechanism (co2lab's postprocessStates). The CO2 and brine
# properties are JutulDarcy's, not co2lab's.
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/johansen_ve.jl
#
# Input: the Johansen data set's sector model, read directly — NPD5.grdecl
# (ECLIPSE grid) and NPD5_Permeability.txt, NPD5_Porosity.txt, in
# johansen_data_directory() (../mrst-data/Johansen; from MRST's Johansen data
# set, https://www.sintef.no/contentassets/124f261f170947a6bc51dd76aea66129/
# johansen.zip). No MRST output is needed.
#
# Plotting backend: GLMakie (a window per figure, 3D views rotatable with the
# mouse; close it to continue) by default; JUTULCO2LAB_BACKEND=CairoMakie
# writes the figures to files instead.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Jutul, JutulDarcy, Printf

year = si_unit(:year)

# ------------------------------ 1. The formation ------------------------------

# co2lab's Johansen VE grid (makeJohansenVEgrid), built from the deck: the
# sector model's cells with permeability above 0.11 mD, columns down to the
# first hole (topSurfaceGrid's discard_below_holes), the rock averaged over
# each column. Across the faults where the caprock breaks off (146 faces) the
# flux compares the two sides' pressures at a common depth (a virtual cell);
# co2lab joins the caprocks smoothly instead (`fault_flux = :co2lab`).
# The formation at the example's temperature, 94 °C, for the fluid tables.
johansen = johansen_ve_domain(; temperature = 273.15 + 94)
mesh, columns, ve, sides = johansen.mesh, johansen.columns, johansen.ve, johansen.sides
@printf("%d columns of %.0f-%.0f m\n", length(ve.geometry.H),
        extrema(ve.geometry.H)...)

# --------------------------- 2. A capillary fringe ----------------------------

# CO2 and brine from JutulDarcy's CO2–brine tables, in pressure and the
# columns' temperature (co2lab's example linearises its own equation of state
# at 300 bar instead); rock compressibility 4.35e-5 /bar at 300 bar.
srw, src, pe = 0.27, 0.2, 5e3

# Below the caprock the capillary pressure falls off linearly with depth, and
# the rock's capillary curve turns it into a CO2 saturation that fades from
# 1 − srw to 0 over a fringe. The curve: entry pressure 5 kPa, brine down to
# its residual 0.27; CO2 relperm zero below the residual CO2 saturation 0.2.
fringe = CapillaryFringeTable(; invPc3D = p -> (1 - srw) * (pe / max(p, pe))^2 + srw,
                              kr3D = s -> max((s - src) / (1 - src), 0.0)^2)
fluid = co2lab_ve_fluid(; vertical = fringe, residual = (srw, src))
rock_c = (reference_pressure = 300e5, compressibility = 4.35e-5 / 1e5)

# ------------------------- 3. Boundary, well, schedule -------------------------

# Parts of the outer edges are open, as co2lab's example sets them
# (makeJohansenVEgrid's searchForBoundaryFaces calls, by logical i and j;
# MRST's 'BACK' is Jutul's :lower, 'FRONT' :upper), the rest closed; the
# open sides at hydrostatic pressure: that of the tables' brine at rest
# (lighter than co2lab's 1000 kg/m³ at 94 °C).
hydrostatic = ve_hydrostatic_profile(fluid, ve)
open_on(dir; i = nothing, j = nothing) =
    ve_select_sides(sides, columns; directions = (dir,), i, j)
open_sides = sort(union(open_on(:lower; i = 1:94, j = 1:4), open_on(:left; i = 1:20),
                        open_on(:right; j = 90:100),
                        open_on(:upper; i = 1:42, j = 50:100)))
bc = ve_boundary_conditions(ve, sides, open_sides; pressure = hydrostatic)

# One injector far down-dip near the major sealing fault, in the column of
# pillar (48, 48), 3.5 Mt of CO2 a year for 50 years; then 950 years of
# migration in 10-year report steps.
well_column = findfirst(==((48, 48)), columns.pillar)
well = setup_ve_well(ve, well_column; name = :inj)
# (the rate as a volume at the fluid's reference CO2 density)
inject = co2_injector_control(TotalRateTarget(3.5e9 / year / fluid.rho_ref[2]), fluid)
model = setup_co2lab_model(ve.domain; wells = [well], fluid,
                           rock_compressibility = rock_c)
f_inj = setup_reservoir_forces(model; control = Dict(:inj => inject), bc = bc)
f_shut = setup_reservoir_forces(model; control = Dict(:inj => DisabledControl()),
                                bc = bc)
dt = [fill(year, 50); fill(10year, 95)]
forces = [k <= 50 ? f_inj : f_shut for k in eachindex(dt)]

# ---------------------------------- 4. Run ------------------------------------

p0 = hydrostatic.(ve.geometry.z_top)
state0 = setup_co2lab_state(model; pressure = p0, forces = f_inj)
prm = setup_parameters(model)
run = JutulCase(model, dt, forces; state0 = state0, parameters = prm)
sim, cfg = setup_co2lab_simulator(run; info_level = -1)
t = @elapsed result = simulate_reservoir(run; simulator = sim, config = cfg)
@printf("1000 years simulated in %.0f s\n", t)

# --------------------------- 5. The plume in 3D -------------------------------

# The caprock with its traps, and over it the columns holding CO2 at the end
# of injection and after 1000 years, coloured by the CO2 saturation's
# capillary height (the depth of the fringe's foot below the caprock).
s = trap_surface(mesh, columns)
traps = trap_analysis(s; method = :node)
for (k, label) in ((50, "end of injection (50 years)"), (length(dt), "1000 years"))
    h = result.states[k][:PlumeHeights][1, :]
    plume = findall(>(0.01), h)
    fig = Figure(size = (1200, 800))
    ax = ve_axis3(fig[1, 1], s; view = (-63, 50), exaggeration = 5, title = label)
    plot_traps!(ax, s, traps; rivers = false)
    plot_ve_surface!(ax, (; s..., polygons = s.polygons[plume]); color = h[plume],
                     colormap = :plasma, offset = -3.0, shading = true, strokewidth = 0)
    plot_ve_map!(ax, s; wells = [well_column], outline = false)
    Colorbar(fig[1, 2]; limits = extrema(h[plume]), colormap = :plasma,
             label = "capillary height (m)")
    show_figure(fig, "johansen_ve_$(k)")
end

# ------------------------- 6. The trapping inventory --------------------------

# Where the CO2 is, by trapping mechanism (co2lab's postprocessStates): in
# structural traps (from the spill-point analysis), residually trapped, still
# free, or out through the open boundary. In the fringe, the split between
# zones follows the fringe's saturation profile.
history = trapping_history(result, model, fluid, s, traps)
@printf("\nafter 1000 years (Mt):\n")
for k in TRAPPING_CATEGORIES
    @printf("  %-20s %7.2f\n", k, history.inventory[end][k] / 1e9)
end
@printf("  %-20s %7.2f\n", :exited, history.exited[end] / 1e9)
fig = plot_trapping_history(history; title = "Johansen: trapping inventory")
show_figure(fig, "johansen_ve_inventory")
