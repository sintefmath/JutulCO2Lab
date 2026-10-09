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

On a 3D axis (`Axis3`, e.g. [`ve_axis3`](@ref)) `plot_ve_map!` draws on the
caprock at its depths, lit: `background` on the surface, depth contours at
their depths, `regions` as an overlay, wells as vertical lines; `plumes`
only on map views.

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
On a 3D axis (`Axis3`) `plot_atlas_map!` draws the map as a surface with its
values as height (co2lab's `surf`). `plot_atlas_map` returns a new `Figure`,
`plot_atlas_map!` the axis.
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
                     colormap = :viridis, shading = false, strokecolor = :black,
                     strokewidth = 1, smooth = true) -> ax

The top surface of a VE domain ([`ve_top_surface`](@ref)) in 3D: the
columns' top faces at their node depths (plus `offset`, e.g. to lift it off
the 3D grid), with their edges (none for `strokewidth = 0`). `color` is one
colour, one value per column (drawn in `colormap`) or one colour per
column; `shading` lights the surface, `smooth` across column edges (normals
shared at the nodes; `false`: each column lit as a facet of its own; colours
stay per column either way). Modifies: `ax` (an `Axis3`). Returns `ax`.
"""
function plot_ve_surface! end

"""
    plot_traps(s, res; title = "", colorbar = true, kwarg...) -> Figure
    plot_traps!(ax, s, res; regions = true, trap_values = nothing,
                colormap = :viridis, rivers = true, labels = false,
                outline = true, faults = true) -> ax

Map view of a trapping analysis — the trap surface `s`
([`trap_surface`](@ref)) and [`trap_analysis`](@ref)'s result `res`, as the
static plots of co2lab's spill-point examples:

- `regions`: each trap's spill region in a light colour (the region
  draining out of the domain in grey);
- the traps in their region's full colour, or coloured by `trap_values`
  (one value per trap, e.g. volumes) in `colormap`;
- `rivers`: the spill paths through the column centroids;
- `labels`: trap numbers at the traps' top columns;
- `outline`: the domain's lateral boundary;
- `faults`: the sealing faults (`s.sealed`) as thick red lines.

On a 3D axis (`Axis3`, e.g. [`ve_axis3`](@ref)) `plot_traps!` colours the
caprock at its depths, lit, with rivers, labels, outline and faults on it.

`plot_traps` returns a new `Figure` (a colour bar when `trap_values` is
given); `plot_traps!` draws into the axis `ax` and returns it. Modifies
(`!` version): `ax`. `s` and `res` are only read.
"""
function plot_traps end
function plot_traps! end

"""
    plot_trapping_history(history; title = "", kwarg...) -> Figure
    plot_trapping_history!(ax, history; legend = true) -> ax

The trapping inventory over a run ([`trapping_history`](@ref)'s result) as
stacked areas, as co2lab's `plotTrappingDistribution`: time in years, mass
in Mt, co2lab's categories, order and colours, the CO2 that has left the
domain on top; categories that stay zero (dissolution, subscale) are left
out of the legend. `plot_trapping_history` returns a new `Figure`,
`plot_trapping_history!` draws into `ax` and returns it. Modifies (`!`
version): `ax`. `history` is only read.
"""
function plot_trapping_history end
function plot_trapping_history! end

"""
    interactive_traps(s; res = nothing, method = :node, traps = true,
                      spill_regions = false, lighting = false, contours = false,
                      path_colours = true, injection_column = 0, exaggeration = 10,
                      view = (-37.5, 30.0), title = "", size = (1700, 1000),
                      widgets = extra -> nothing, charts = true, selection = true,
                      on_click = nothing, on_draw = nothing) -> NamedTuple
    interactive_traps!(viewer, s; res = nothing, title = "", injection_column = 0)
        -> viewer

The interactive trap viewer, co2lab's `interactiveTrapping`: the trap
surface `s` ([`trap_surface`](@ref)) in 3D with its traps, rivers, outline
and sealing faults, and a panel of widgets over the window's right side.
Interactive with GLMakie (a window) or WGLMakie (a browser page); other
backends draw the first view only.

- A left click on the caprock selects a column (co2lab: the injection
  point; marked yellow, with a red vertical line): the trap it drains into
  and every trap that trap spills into are coloured in spill order, bright
  red to nearly black (grey with the "Red path" switch off), their outline
  blue and their rivers thick; with shift (or the middle button) the traps
  spilling into it instead (outline green). A right click clears the
  selection.
- The panel: the trap analysis's method (recomputed when changed), the
  vertical exaggeration; switches for the traps, the spill regions (off: the
  caprock coloured by depth), lighting, ten depth contours and the red path;
  buttons for a view from above, the starting view and the selected trap
  alone (its floor, filled to the spill depth, red; its ceiling blue; the
  caprock around faint; again to return); the selection's numbers; co2lab's
  pie (the trap's volume, the volume further along the path, the other
  traps' volume) and its bars (each path trap's volume, log scale).
- `res`: a [`trap_analysis`](@ref) result of `s` to show (default: computed
  with `method`). The keywords set the widgets' starting state as co2lab's
  options do (`path_colours` is its `colorpath`, `injection_column` its
  `injpt`); `view` is the starting (azimuth, elevation) in MATLAB's degrees,
  as [`ve_axis3`](@ref).
- `widgets`: a function given an empty `GridLayout` at the panel's top, to
  make a caller's own widgets in (e.g. a formation menu); called before the
  viewer makes its own, so that a caller's dropdown list takes the clicks
  over the widgets below it (Makie gives a click to the first-made of
  overlapping widgets). Its return value is the viewer's `widgets`.
- `charts = false`: no pie and bars (room for a caller's widgets);
  `selection = false`: clicks select nothing, and the "Red path" switch and
  the "Trap detail" button are left out (`toggles` without `red`, `detail`
  `nothing`).
- `on_click(column, button)`: called for each click on the caprock (the
  column under the pointer, 0 if none; `Mouse.left`, `Mouse.right`, ...),
  instead of selecting traps (`selection` off then).
- `on_draw(viewer)`: called after each drawing of the caprock (a new
  surface, the lighting switched), which removes all plots of the axis: for
  a caller's own plots on it.
- [`show_column_values!`](@ref) colours the caprock by values per column
  instead (e.g. capacities), with a colour bar.
- Differences from co2lab: the method asked for is used (co2lab's
  `interactiveTrapping` always uses the node-based one; `co2lab-issues.org`);
  the default method is `:node`, what co2lab shows; the trap figure is shown
  in the same axis; the bars are in spill order and coloured as the map; no
  "Setup VE simulation" button yet (roadmap Stage 11 item 3).

`interactive_traps` returns a new NamedTuple, the viewer: `fig` (show it
with `display` or [`show_figure`](@ref)), `ax` (the `Axis3`), `title`,
`panel`, `backdrop`, `extra` (the `GridLayout` given to `widgets`),
`widgets` (what that function returned), the widgets `method`,
`exaggeration`, `toggles` (`traps`, `regions`, `light`, `contours`, `red`),
`top`, `reset`, `detail`, `info`, the charts' areas `pie` and `bars` and
their plots `charts` (empty without charts), the `colorbar`, `on_click`,
`on_draw`, and its state: `values` (a `Ref` to the values shown, or
`nothing`), `selectable`, `view`, `current` (a `Ref`
to `(s, res, vol, title)`), `selection` (a `Ref` to `(column, upstream)`),
`plots`, `in_detail`, `lock`, `pressed_at`. `interactive_traps!` shows the surface
`s` in an existing viewer (e.g. another formation), with the viewer's
method; it returns `viewer`. Modifies (`!` version): `viewer` (its figure's
plots and its state). `s` and `res` are only read.
"""
function interactive_traps end
function interactive_traps! end

"""
    show_column_values!(viewer, values; label = "", colormap = :viridis) -> viewer

Colours the caprock of an [`interactive_traps`](@ref) viewer by `values`,
one per column (e.g. capacities, depths), in `colormap` over their finite
range (`NaN`: light grey), with a colour bar labelled `label` at the bottom
left; the traps' and spill regions' colours give way to them, a selection
keeps its outline, rivers and marker. `values = nothing` returns to the
viewer's own colours; a new surface (`interactive_traps!`) does too.
Modifies: `viewer` (its colours, colour bar and `values`). Returns `viewer`.
`values` is only read (copied).
"""
function show_column_values! end

"""
    ve_axis3(pos, surface; exaggeration = 10.0, view = (-37.5, 30.0), title = "")
        -> Axis3

A 3D axis for a VE top surface (e.g. [`trap_surface`](@ref),
[`ve_top_surface`](@ref)) at figure position `pos` (e.g. `fig[1, 1]`): depth
downwards, the surface's extent at equal lateral scale and its depth range
`exaggeration` times enlarged, seen from `view = (az, el)` in degrees as
MATLAB's `view` (so co2lab's views carry over). Rotatable with the mouse in
GLMakie. To set its limits, give the depth bounds deepest first
(`limits!(ax, x1, x2, y1, y2, z_deep, z_shallow)`): Makie's `limits!` with
the shallow bound first turns the depth axis upwards. Modifies: the figure of
`pos` (the axis added). Returns the new `Axis3`. `surface` is only read.
"""
function ve_axis3 end

"""
    show_figure(fig, name; dir = "figures") -> nothing

Show `fig` in a window of its own with GLMakie — in a script (not the REPL)
waiting until the window is closed, since Julia closes its windows on exit —
or save it as `dir/name.png` with CairoMakie, so that example scripts run
either way. Other backends: `display(fig)`. Returns `nothing`.
"""
function show_figure end
