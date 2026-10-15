# co2lab-spillpoint's trapsJohansen: structural trapping in the Johansen
# formation, proposed for CO2 from Mongstad and used in the Stuttgart
# benchmark study (Eigestad et al.), estimated on three models — the NPD's
# sector model (the benchmark's), the NPD's full-field model and the CO2
# Storage Atlas — from the benchmark's injection point.
#
#   julia --project=examples examples/trapping/traps_johansen.jl
#
# Needs the NPD decks in johansen_data_directory() (downloadable from SINTEF:
# the error message names the files). Plotting backend: GLMakie (a window per
# figure, 3D views rotatable with the mouse; close it to continue) by default;
# JUTULCO2LAB_BACKEND=CairoMakie writes the figures to files instead. The
# traps reached from the injection point are drawn on the caprock, then each
# model is shown in the interactive trap viewer, as in co2lab.

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

# The Atlas's caprock coloured by depth, the NPD models' outlines at their
# depths, and the well.
atlas = last(surfaces[3])
fig = Figure(size = (1200, 900))
ax = ve_axis3(fig[1, 1], atlas; view = (-62, 60), exaggeration = 5,
              title = "Johansen: Atlas caprock depth, NPD models outlined")
plot_ve_map!(ax, atlas; background = atlas.zc, outline = false)
plot_ve_outline!(ax, last(surfaces[2]); color = :red, label = "NPD full field")
plot_ve_outline!(ax, last(surfaces[1]); color = :black, label = "NPD sector")
zmin, zmax = extrema(atlas.z)
lines!(ax, [Point3f(w..., zmin), Point3f(w..., zmax)]; color = :blue, linewidth = 3,
       label = "injection well")
Legend(fig[2, 1], ax; orientation = :horizontal)
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
    a = ve_axis3(fig[1, k], s; view = (-80, 64), exaggeration = 5,
                 title = "$name: traps reached (red)")
    plot_traps!(a, s, res; rivers = false)
    plot_ve_map!(a, s; regions = res.traps .∈ Ref(reached), region_color = :red,
                 wells = [col], outline = false)
end
show_figure(fig, "johansen_traps_from_well")

# ----------------------- 3. The interactive trap viewer -----------------------

# Each model in the interactive viewer (interactive_traps, co2lab's
# interactiveTrapping), as co2lab opens them: cell-based traps, lit, spill
# regions shown, the path traps in grey, the well's column selected. Close a
# window to go on to the next model.
for (name, s) in surfaces
    v = interactive_traps(s; method = :cell, lighting = true, spill_regions = true,
                          path_colours = false, injection_column = nearest(s),
                          view = (-80, 64), title = "Johansen, $name")
    show_figure(v.fig, "johansen_interactive_$(replace(name, ' ' => '_'))")
end
