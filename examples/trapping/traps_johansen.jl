# co2lab-spillpoint's trapsJohansen: structural trapping in the Johansen
# formation, proposed for CO2 from Mongstad and used in the Stuttgart
# benchmark study (Eigestad et al.), estimated on three models — the NPD's
# sector model (the benchmark's), the NPD's full-field model and the CO2
# Storage Atlas — from the benchmark's injection point.
#
#   scripts/julia --project=JutulCO2Lab/examples JutulCO2Lab/examples/trapping/traps_johansen.jl
#
# Needs the NPD decks in johansen_data_directory() (downloadable from SINTEF:
# the error message names the files). Plotting backend: GLMakie (a window per
# figure; close it to continue) by default; JUTULCO2LAB_BACKEND=CairoMakie
# writes the figures to files instead. co2lab shows each model in its
# interactive trap viewer (roadmap Stage 11); here the traps reached from the
# injection point are drawn on a map.

using JutulCO2Lab
backend = Symbol(get(ENV, "JUTULCO2LAB_BACKEND", "GLMakie"))
@eval using $backend
using Printf

# ---------------------------- 1. The three models -----------------------------

# The NPD models are faulted corner-point grids. Their faults are sealing by
# default: closed boundaries of the trap analysis, so CO2 collects against
# them (drawn red on the maps). co2lab instead stitches the top surface
# across faults, as if they were open (`faults = :open`).
models = ["NPD sector" => johansen_mesh(:sector),
          "NPD full field" => johansen_mesh(:full_field),
          "Atlas" => atlas_mesh("Johansenfm")]
surfaces = [name => trap_surface(m.mesh, m.columns) for (name, m) in models]

# The trapping capacity with sealing and with open faults.
for (name, m) in models[1:2]
    capacity(f) = (s = trap_surface(m.mesh, m.columns; faults = f);
                   sum(trap_volumes(s, trap_analysis(s))))
    @printf("%-15s capacity with sealing faults %.2e m3, with open faults %.2e m3\n",
            name, capacity(:sealing), capacity(:open))
end

# The injection well of the benchmark, and on each model the column nearest
# to it.
w = johansen_well_point()
nearest(s) = argmin((s.cx .- w[1]) .^ 2 .+ (s.cy .- w[2]) .^ 2)

# The three models' outlines on the Atlas's caprock depth, and the well.
atlas = last(surfaces[3])
fig = plot_ve_map(atlas; background = atlas.zc, depth_contours = 0,
                  title = "Johansen: Atlas caprock depth (m), NPD models outlined")
ax = current_axis()
plot_ve_outline!(ax, last(surfaces[2]); color = :red, label = "NPD full field")
plot_ve_outline!(ax, last(surfaces[1]); color = :black, label = "NPD sector")
scatter!(ax, [w[1]], [w[2]]; color = :blue, markersize = 14, label = "injection well")
axislegend(ax; position = :rb)
show_figure(fig, "johansen_models")

# ----------------- 2. Traps reached from the injection point ------------------

# Cell-based traps on each model. CO2 injected at the well rises to the trap
# whose spill region holds the well's column, fills it and spills onwards:
# the traps it can fill are that trap and the ones downstream of it.
fig = Figure(size = (1800, 750))
for (k, (name, s)) in enumerate(surfaces)
    res = trap_analysis(s; method = :cell)
    volumes = trap_volumes(s, res)
    col = nearest(s)
    first_trap = res.trap_regions[col]
    reached = first_trap == 0 ? Int[] : downstream_traps(res, first_trap)
    @printf("%-15s %6d columns, %4d traps, capacity %.2e m3; ", name,
            length(s.zc), length(volumes), sum(volumes))
    @printf("from the well %d traps, %.2e m3\n", length(reached),
            sum(volumes[reached]; init = 0.0))

    # all traps and their spill regions; the traps reached in red, the well
    a = Axis(fig[1, k]; title = "$name: traps reached (red)", aspect = DataAspect(),
             xticks = LinearTicks(3))
    plot_traps!(a, s, res; rivers = false)
    plot_ve_map!(a, s; regions = res.traps .∈ Ref(reached), region_color = :red,
                 depth_contours = 0, wells = [col])
end
show_figure(fig, "johansen_traps_from_well")
