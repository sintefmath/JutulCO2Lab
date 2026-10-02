# co2lab's exampleVEJohansen.m with JutulCO2Lab: CO2 injected into the
# Johansen formation (the NPD sector model) for 50 years and followed for 950
# more, in a vertical-equilibrium model with a capillary fringe — the CO2 and
# brine meet in a transition zone set by the rock's capillary pressure, not at
# a sharp interface — and residual trapping. At the end, where the CO2 is, by
# trapping mechanism (co2lab's postprocessStates).
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/johansen_ve.jl
#
# Input: co2lab's grid and open boundary (analysis/reference-data/
# johansen_grids.mat and stage1_johansen.mat, from scripts/mrst_dump_stage1.m)
# and the example's CO2 properties at its reference conditions
# (johansen_case.mat, from scripts/mrst_ref_johansen.m 'case').
#
# Plotting backend: GLMakie (a window per figure, 3D views rotatable with the
# mouse; close it to continue) by default; JUTULCO2LAB_BACKEND=CairoMakie
# writes the figures to files instead.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Jutul, JutulDarcy, MAT, Printf

refdir = joinpath(@__DIR__, "..", "..", "analysis", "reference-data")
year = 365.25 * 86400

# ------------------------------ 1. The formation ------------------------------

# co2lab's Johansen VE grid: the sector model's cells with permeability above
# 0.11 mD, columns down to the first hole (topSurfaceGrid's
# discard_below_holes), the rock averaged over each column.
d = matread(joinpath(refdir, "johansen_grids.mat"))
mesh = ve_mesh_from_mrst(d["G3"])
columns, _ = ve_columns(mesh)
K = d["rock3D"]["perm"]
ve = ve_domain(mesh, columns; perm_x = K[:, 1], perm_y = K[:, 2],
               poro = vec(d["rock3D"]["poro"]))
@printf("%d columns of %.0f-%.0f m\n", length(ve.geometry.H),
        extrema(ve.geometry.H)...)

# --------------------------- 2. A capillary fringe ----------------------------

# The example's constants: brine 1000 kg/m³, CO2 at 300 bar and 94 °C from
# co2lab's equation of state, linearised (compressible CO2, incompressible
# brine), rock compressibility 4.35e-5 /bar.
c = matread(joinpath(refdir, "johansen_case.mat"))
srw, src, pe = c["srw"], c["src"], c["pe"]

# Below the caprock the capillary pressure falls off linearly with depth, and
# the rock's capillary curve turns it into a CO2 saturation that fades from
# 1 − srw to 0 over a fringe. The curve: entry pressure 5 kPa, brine down to
# its residual 0.27; CO2 relperm zero below the residual CO2 saturation 0.2.
fringe = CapillaryFringeTable(; invPc3D = p -> (1 - srw) * (pe / max(p, pe))^2 + srw,
                              kr3D = s -> max((s - src) / (1 - src), 0.0)^2)
fluid = co2lab_ve_fluid(; default_pvt = nothing,
                        vertical = fringe, residual = (srw, src),
                        co2_rho_ref = c["co2_rho"], wat_rho_ref = c["rhow"],
                        co2_rho_pvt = (c["co2_c"], c["p_ref"]),
                        wat_rho_pvt = (c["wat_c"], c["p_ref"]),
                        co2_mu_ref = c["muco2"], wat_mu_ref = c["muw"])
rock_c = (reference_pressure = c["p_ref"], compressibility = c["c_rock"])

# ------------------------- 3. Boundary, well, schedule -------------------------

# Parts of the outer edges are open, as co2lab's example sets them
# (makeJohansenVEgrid's bcIxVE: its faces, matched to our sides by column and
# edge midpoint), the rest closed; the open sides at hydrostatic pressure.
g, rhow = c["gravity"], c["rhow"]
s1 = matread(joinpath(refdir, "stage1_johansen.mat"))
sides, _ = ve_boundary_sides(mesh, columns, ve.domain)
Nt = Int.(d["Gt"]["faces"]["neighbors"])
ct = d["Gt"]["faces"]["centroids"]
open_faces = Set(Int.(vec(s1["bc_faces_ve"])))
function is_open(i)
    col, m = sides.column[i], sides.edge_midpoint[i]
    return any(f in open_faces && 0 in Nt[f, :] && col in Nt[f, :] &&
               hypot(ct[f, 1] - m[1], ct[f, 2] - m[2]) < 1e-6 for f in axes(Nt, 1))
end
bc = ve_boundary_conditions(ve, sides, filter(is_open, 1:length(sides));
                            pressure = z -> rhow * g * z)

# One injector far down-dip near the major sealing fault, 3.5 Mt of CO2 a
# year for 50 years; then 950 years of migration in 10-year report steps.
well_column = Int(only(s1["well_cells"]))
well = setup_ve_well(ve, well_column; name = :inj)
inject = co2_injector_control(TotalRateTarget(c["inj_rate"]), fluid)
model = setup_co2lab_model(ve.domain; wells = [well], fluid,
                           rock_compressibility = rock_c)
f_inj = setup_reservoir_forces(model; control = Dict(:inj => inject), bc = bc)
f_shut = setup_reservoir_forces(model; control = Dict(:inj => DisabledControl()),
                                bc = bc)
dt = [fill(year, 50); fill(10year, 95)]
forces = [k <= 50 ? f_inj : f_shut for k in eachindex(dt)]

# ---------------------------------- 4. Run ------------------------------------

p0 = ve_hydrostatic_pressure(ve.geometry; rho_w = rhow, g = g)
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
