# Plotting functions, implemented in the Makie extension
# (ext/JutulCO2LabMakieExt): available once a Makie backend is loaded —
# `using GLMakie` (interactive windows) or `using CairoMakie` (files, no
# display needed). 3D views use Jutul's own `plot_cell_data` and `plot_mesh`.

"""
    plot_ve_map(surface; kwarg...) -> Figure
    plot_ve_map!(ax, surface; background = nothing, colormap = :viridis,
                 depth_contours = 40, plumes = [], plume_threshold = 0.3,
                 plume_colors = [:red, :blue, :cyan, :green, :magenta],
                 labels = nothing, wells = Int[], regions = nothing,
                 region_color = (:red, 0.3), outline = true, colorbar = true)

Map view of a VE domain, as co2lab's `mapPlot`, on the exact column polygons
of `surface` ([`ve_top_surface`](@ref)) rather than a resampled raster:

- `background`: a field per column, drawn as coloured polygons;
- `depth_contours`: number of caprock depth contour lines (0 for none),
  traced on the surface's node depths;
- `plumes`: fields per column (e.g. plume heights), each outlined where it
  crosses `plume_threshold`, in `plume_colors`, with `labels` as legend;
- `wells`: columns marked as wells;
- `regions`: a per-column mask (e.g. traps) filled in `region_color`;
- `outline`: the domain's lateral boundary.

`plot_ve_map` returns a new `Figure`; `plot_ve_map!` draws into the axis
`ax` and returns it. Modifies (`!` version): `ax`. `surface` is only read.
"""
function plot_ve_map end
function plot_ve_map! end

"""
    plot_atlas_map(map; kwarg...) -> Figure
    plot_atlas_map!(ax, map; style = :heatmap, levels = 10, colormap = :viridis)

An Atlas map as read by [`read_atlas_maps`](@ref) (`(name, variant, meta,
data)`), at the positions of its values: `style = :heatmap` colours the
values, `:contour` draws `levels` contour lines (co2lab's `contourAtlas`).
`plot_atlas_map` returns a new `Figure`, `plot_atlas_map!` the axis.
Modifies (`!` version): `ax`. `map` is only read.
"""
function plot_atlas_map end
function plot_atlas_map! end

"""
    plot_atlas_context!(ax; coastline = true, wells = true) -> ax

Norway's coastline ([`atlas_coastline`](@ref)) and the shelf's wells
([`atlas_wells`](@ref)) for orientation on a map axis. Modifies: `ax`.
Returns `ax`.
"""
function plot_atlas_context! end

"""
    plot_ve_outline!(ax, surface; color = :black, linewidth = 2, label = nothing)
        -> ax

The lateral boundary of a VE domain (`surface.outline` of
[`ve_top_surface`](@ref)) as lines. Modifies: `ax`. Returns `ax`.
"""
function plot_ve_outline! end

"""
    plot_ve_surface!(ax, surface; color = (:steelblue, 0.6), offset = 0.0,
                     strokecolor = :black, strokewidth = 1) -> ax

The top surface of a VE domain ([`ve_top_surface`](@ref)) in 3D: the
columns' top faces at their node depths (plus `offset`, e.g. to lift it off
the 3D grid), with their edges. Modifies: `ax` (an `Axis3`). Returns `ax`.
"""
function plot_ve_surface! end

"""
    show_figure(fig, name; dir = "figures") -> nothing

Show `fig` in a window of its own with GLMakie — in a script (not the REPL)
waiting until the window is closed, since Julia closes its windows on exit —
or save it as `dir/name.png` with CairoMakie, so that example scripts run
either way. Other backends: `display(fig)`. Returns `nothing`.
"""
function show_figure end
