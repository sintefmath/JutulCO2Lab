# Plotting functions. They are declared and documented here, and implemented in
# the Makie extension (ext/JutulCO2LabMakieExt), which loads with any Makie
# backend: `using GLMakie` (interactive windows) or `using CairoMakie` (image
# files, no display needed). 3D views of fine grids use Jutul's own
# `plot_cell_data` and `plot_mesh`.
#
# Each plot comes in two forms: `plot_x(...)` makes a new `Figure`, and
# `plot_x!(ax, ...)` draws into an existing axis. A 2D `Axis` gives a map
# view; an `Axis3` (e.g. from `ve_axis3`) draws on the caprock in 3D.

# ---------------------------------- Map views ---------------------------------

"""
    plot_ve_map(surface; title = "", colorbar = true, kwarg...) -> Figure
    plot_ve_map!(ax, surface; background = nothing, colormap = :viridis,
                 depth_contours = 40, plumes = [], plume_threshold = 0.3,
                 plume_colors = [:red, :blue, :cyan, :green, :magenta],
                 labels = nothing, wells = Int[], regions = nothing,
                 region_color = (:red, 0.3), outline = true) -> ax

Map view of a VE domain, co2lab's `mapPlot`. It is drawn on the exact column
polygons of the top surface, not on a resampled raster.

# Arguments

- `ax`: the axis to draw into: a 2D `Axis` (map view) or an `Axis3` (on the
  caprock in 3D, see Notes).
- `surface`: the domain's top surface ([`ve_top_surface`](@ref)).
- `background`: one value per column, drawn as coloured polygons in
  `colormap`; `nothing` for none.
- `depth_contours`: number of caprock depth contour lines, traced on the
  surface's node depths; 0 for none. Default 40 on a map, 0 on an `Axis3`.
- `plumes`: fields with one value per column (e.g. plume heights); each is
  outlined where it crosses `plume_threshold`, in `plume_colors` (at most
  one field per colour).
- `labels`: legend labels, one per plume; `nothing` for no legend.
- `wells`: columns to mark as wells.
- `regions`: a per-column mask (e.g. traps), filled in `region_color`.
- `outline`: draw the domain's lateral boundary.
- `title`, `colorbar` (`plot_ve_map` only): the axis title, and a colour bar
  for `background`. Other keywords are passed on to `plot_ve_map!`.

# Returns

`plot_ve_map`: a new `Figure` with one map axis. `plot_ve_map!`: `ax`.

Modifies: `ax` (`plot_ve_map!`)

# Notes

- On an `Axis3` the map lies on the caprock at its depths, lit: `background`
  on the surface, contours at their depths, `regions` as an overlay, wells
  as vertical lines. `plumes` are only drawn on map views (an error
  otherwise).
"""
function plot_ve_map end

"""
    plot_ve_map!(ax, ...)

The in-place form of [`plot_ve_map`](@ref), documented there.
"""
function plot_ve_map! end

"""
    plot_atlas_map(map; title = "<name> <variant>", kwarg...) -> Figure
    plot_atlas_map!(ax, map; style = :heatmap, levels = 10, colormap = :viridis)
        -> ax

A raw CO2 Storage Atlas map, drawn at the positions of its values, as
co2lab's `contourAtlas`.

# Arguments

- `ax`: the axis to draw into: a 2D `Axis`, or an `Axis3` for a surface with
  the map's values as height (co2lab's `surf`; `style`, `levels` ignored).
- `map`: one map from [`read_atlas_maps`](@ref), a `(name, variant, meta,
  data)` NamedTuple.
- `style`: `:heatmap` colours the values; `:contour` draws contour lines.
- `levels`: for `:contour`, the number of evenly spaced contour levels, or
  the levels themselves.
- `colormap`: colour map of the heatmap or contour lines.
- `title` (`plot_atlas_map` only): the axis title; default the map's name
  and variant. Other keywords are passed on to `plot_atlas_map!`.

# Returns

`plot_atlas_map`: a new `Figure` (with a colour bar for a heatmap).
`plot_atlas_map!`: `ax`.

Modifies: `ax` (`plot_atlas_map!`)
"""
function plot_atlas_map end

"""
    plot_atlas_map!(ax, ...)

The in-place form of [`plot_atlas_map`](@ref), documented there.
"""
function plot_atlas_map! end

"""
    plot_atlas_context!(ax; coastline = true, wells = true) -> ax

Adds Norway's coastline ([`atlas_coastline`](@ref)) and the shelf's wells
([`atlas_wells`](@ref)) to a map axis, for orientation. Returns `ax`.

Modifies: `ax`
"""
function plot_atlas_context! end

"""
    plot_ve_outline!(ax, surface; color = :black, linewidth = 2, label = nothing)
        -> ax

Draws the lateral boundary of a VE domain (`surface.outline` of
[`ve_top_surface`](@ref)) as lines; on an `Axis3`, just above the caprock.
`label` is its legend entry. Returns `ax`.

Modifies: `ax`
"""
function plot_ve_outline! end

# ---------------------------------- 3D surface --------------------------------

"""
    plot_ve_surface!(ax, surface; color = (:steelblue, 0.6), offset = 0.0,
                     colormap = :viridis, shading = false, strokecolor = :black,
                     strokewidth = 1, smooth = true) -> ax

Draws the top surface of a VE domain in 3D: each column's top face at its
node depths, with its edges, as co2lab draws its top-surface grid.

# Arguments

- `ax`: an `Axis3` (e.g. [`ve_axis3`](@ref)).
- `surface`: the top surface ([`ve_top_surface`](@ref),
  [`trap_surface`](@ref)).
- `color`: one colour for all; one value per column, drawn in `colormap`;
  or one colour per column.
- `offset`: added to the depths (m), e.g. to lift the surface off a 3D grid.
- `shading`: light the surface.
- `smooth`: blend the lighting across column edges (normals shared at the
  nodes); `false` lights each column as a facet of its own. Colours stay
  per column either way.
- `strokecolor`, `strokewidth`: the column edges; none for `strokewidth = 0`.

# Returns

`ax`.

Modifies: `ax`
"""
function plot_ve_surface! end

# ----------------------------- Traps and inventory ----------------------------

"""
    plot_traps(s, res; title = "", colorbar = true, kwarg...) -> Figure
    plot_traps!(ax, s, res; regions = true, trap_values = nothing,
                colormap = :viridis, rivers = true, labels = false,
                outline = true, faults = true) -> ax

Map view of a trapping analysis, as the static plots of co2lab's spill-point
examples: spill regions, traps, rivers and sealing faults.

# Arguments

- `ax`: the axis to draw into: a 2D `Axis` (map view) or an `Axis3` (the
  caprock coloured at its depths, lit, with the lines on it).
- `s`: the trap surface ([`trap_surface`](@ref)).
- `res`: its trap analysis ([`trap_analysis`](@ref)).
- `regions`: colour each trap's spill region in a light colour, and the
  region draining out of the domain grey.
- `trap_values`: one value per trap (e.g. volumes) to colour the traps by,
  in `colormap`; `nothing` colours them in their region's full colour.
- `rivers`: the spill paths through the column centroids.
- `labels`: trap numbers at the traps' top columns.
- `outline`: the domain's lateral boundary.
- `faults`: the sealing faults (`s.sealed`) as thick red lines.
- `title`, `colorbar` (`plot_traps` only): the axis title, and a colour bar
  when `trap_values` is given. Other keywords are passed on to `plot_traps!`.

# Returns

`plot_traps`: a new `Figure` with one map axis. `plot_traps!`: `ax`.

Modifies: `ax` (`plot_traps!`)
"""
function plot_traps end

"""
    plot_traps!(ax, ...)

The in-place form of [`plot_traps`](@ref), documented there.
"""
function plot_traps! end

"""
    plot_trapping_history(history; title = "", kwarg...) -> Figure
    plot_trapping_history!(ax, history; legend = true) -> ax

The trapping inventory over a run as stacked areas, as co2lab's
`plotTrappingDistribution`: time in years, mass in Mt, co2lab's categories,
order and colours, with the CO2 that has left the domain on top.

# Arguments

- `ax`: the axis to draw into.
- `history`: the result of [`trapping_history`](@ref).
- `legend`: add a legend; categories that stay zero (e.g. dissolution,
  subscale) are left out of it.
- `title` (`plot_trapping_history` only): the axis title. Other keywords
  are passed on to `plot_trapping_history!`.

# Returns

`plot_trapping_history`: a new `Figure`. `plot_trapping_history!`: `ax`.

Modifies: `ax` (`plot_trapping_history!`)
"""
function plot_trapping_history end

"""
    plot_trapping_history!(ax, ...)

The in-place form of [`plot_trapping_history`](@ref), documented there.
"""
function plot_trapping_history! end

# ------------------------------ Interactive viewer ----------------------------

"""
    interactive_traps(s; res = nothing, method = :node, traps = true,
                      spill_regions = false, lighting = false, contours = false,
                      path_colours = true, injection_column = 0, exaggeration = 10,
                      view = (-37.5, 30.0), title = "", size = (1700, 1000),
                      widgets = extra -> nothing, charts = true, selection = true,
                      on_click = nothing, on_draw = nothing) -> NamedTuple
    interactive_traps!(viewer, s; res = nothing, title = "", injection_column = 0)
        -> viewer

The interactive trap viewer, co2lab's `interactiveTrapping`: the trap surface
in 3D with its traps, rivers, outline and sealing faults, and a panel of
widgets on the right. Interactive with GLMakie (a window) or WGLMakie (a
browser page); other backends draw the first view only.

`interactive_traps` makes a new viewer; `interactive_traps!` shows another
surface (e.g. another formation) in an existing one, with the viewer's
current trap-analysis method.

# Using the viewer

- **Left click** on the caprock selects a column (co2lab's injection point;
  yellow, with a red vertical line). The trap it drains into and every trap
  downstream are coloured in spill order, bright red to nearly black (grey
  with "Red path" off), outlined blue, with thick rivers.
- **Shift-click** (or middle button) shows the traps spilling *into* it
  instead (outlined green). **Right click** clears the selection.
- **Panel**: the trap-analysis method (recomputed on change), vertical
  exaggeration; switches for traps, spill regions (off: caprock coloured by
  depth), lighting, ten depth contours and the red path; buttons for a top
  view, the starting view, and "Trap detail" (the selected trap alone: its
  floor filled to the spill depth red, its ceiling blue, the caprock around
  faint; press again to return); the selection's numbers; co2lab's pie (the
  trap, the traps further along the path, the rest) and bar chart (each path
  trap's volume, log scale).
- [`show_column_values!`](@ref) colours the caprock by any per-column values
  instead (e.g. capacities).

# Arguments

- `s`: the trap surface ([`trap_surface`](@ref)).
- `res`: a [`trap_analysis`](@ref) result of `s` to show; default computed
  with `method` (`:node` or `:cell`).
- `traps`, `spill_regions`, `lighting`, `contours`, `path_colours`: the
  switches' starting states (`path_colours` is co2lab's `colorpath`).
- `injection_column`: the column selected at start (co2lab's `injpt`; 0 for
  none).
- `exaggeration`: starting vertical exaggeration (1, 2, 5, 10, 20, 50 or 100).
- `view`: starting `(azimuth, elevation)` in MATLAB's degrees, as
  [`ve_axis3`](@ref).
- `title`, `size`: the figure's title and size (pixels).
- `widgets`: a function `extra -> w` given an empty `GridLayout` at the
  panel's top, for the caller's own widgets (e.g. a formation menu). It is
  called before the viewer makes its own widgets, so that a caller's
  dropdown list gets the clicks over the widgets below it.
- `charts`: show the pie and bars; `false` leaves room for caller widgets.
- `selection`: `false`: clicks select nothing, and the "Red path" switch and
  "Trap detail" button are left out.
- `on_click`: `(column, button) -> ...`, called for each click on the
  caprock (`column` 0 if none; `button` e.g. `Mouse.left`) instead of
  selecting traps; turns `selection` off.
- `on_draw`: `viewer -> ...`, called after each redraw of the caprock (a new
  surface, lighting switched), which removes all plots of the axis: for the
  caller's own plots on it.

# Returns

`interactive_traps`: the viewer, a new NamedTuple:

- `fig`: the `Figure` (show it with `display` or [`show_figure`](@ref)).
- `ax`: the `Axis3`; `title`: the title `Label`.
- `panel`, `backdrop`: the panel's `GridLayout` and its background `Box`.
- `extra`, `widgets`: the layout given to `widgets`, and what it returned.
- `method`, `exaggeration`: the two `Menu`s.
- `toggles`: NamedTuple of `Toggle`s `traps`, `regions`, `light`,
  `contours`, and `red` (only with selection).
- `top`, `reset`, `detail`: the `Button`s (`detail` is `nothing` without
  selection); `info`: the text `Label`.
- `pie`, `bars`: the chart areas (`nothing` without charts); `charts`: their
  plots (empty without charts).
- `colorbar`: the colour bar of [`show_column_values!`](@ref).
- `on_click`, `on_draw`, `view`: as given; `selectable`: whether clicks
  select traps.
- State: `values` (`Ref` to the values shown, or `nothing`), `current`
  (`Ref` to `(s, res, vol, title)`), `selection` (`Ref` to `(column,
  upstream)`), `plots` (`Dict` of the viewer's plots), `in_detail`
  (`Ref{Bool}`), `lock` (serialises drawing), `pressed_at` (mouse presses).

`interactive_traps!`: `viewer`.

Modifies: `viewer` (`interactive_traps!`: its figure's plots and its state)

# Notes

Differences from co2lab:

- The method asked for is used (co2lab's viewer always uses the node-based
  one); the default `:node` is what co2lab shows.
- The trap figure is shown in the same axis ("Trap detail").
- The bars are in spill order and coloured as the map.
- No "Setup VE simulation" button yet.
"""
function interactive_traps end

"""
    interactive_traps!(viewer, s; ...)

The in-place form of [`interactive_traps`](@ref), documented there.
"""
function interactive_traps! end

"""
    show_column_values!(viewer, values; label = "", colormap = :viridis) -> viewer

Colours the caprock of an [`interactive_traps`](@ref) viewer by one value
per column (e.g. capacities, depths), with a colour bar at the bottom left.

# Arguments

- `viewer`: the viewer.
- `values`: one value per column, drawn in `colormap` over their finite
  range (`NaN`: light grey); copied. `nothing` returns to the viewer's own
  colours, as does showing a new surface ([`interactive_traps!`](@ref)).
- `label`: the colour bar's label.

# Returns

`viewer`.

Modifies: `viewer` (its colours, colour bar and `values`)

# Notes

- The values replace the traps' and spill regions' colours; a selection
  keeps its outline, rivers and marker.
"""
function show_column_values! end

# -------------------------------- Axes, output --------------------------------

"""
    ve_axis3(pos, surface; exaggeration = 10.0, view = (-37.5, 30.0), title = "")
        -> Axis3

A 3D axis suited to a VE top surface: depth pointing down, x and y at equal
scale, depth exaggerated. Rotatable with the mouse in GLMakie.

# Arguments

- `pos`: the figure position to place it at (e.g. `fig[1, 1]`).
- `surface`: the top surface whose extent sets the aspect
  ([`ve_top_surface`](@ref), [`trap_surface`](@ref)).
- `exaggeration`: vertical exaggeration of the depth range.
- `view`: `(azimuth, elevation)` in degrees, as MATLAB's `view`, so that
  co2lab's views carry over. From straight above (elevation ≥ 89) the depth
  axis is hidden.
- `title`: the axis title.

# Returns

A new `Axis3`.

Modifies: the figure of `pos` (the axis is added to it)

# Notes

- To set limits, give the depth bounds deepest first:
  `limits!(ax, x1, x2, y1, y2, z_deep, z_shallow)`. The shallow bound first
  turns the depth axis upwards.
"""
function ve_axis3 end

"""
    show_figure(fig, name; dir = "figures") -> nothing

Shows or saves a figure, so that example scripts run with either backend:

- GLMakie: shows `fig` in a window of its own. In a script (not the REPL)
  it waits until the window is closed, since Julia closes its windows on
  exit.
- CairoMakie: saves it as `dir/name.png` (creating `dir`) and prints the path.
- Other backends: `display(fig)`.

Returns `nothing`.
"""
function show_figure end
