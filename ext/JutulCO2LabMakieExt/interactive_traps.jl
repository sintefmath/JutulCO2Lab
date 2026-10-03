# The interactive trap viewer (roadmap Stage 11 item 5): co2lab's
# interactiveTrapping in Makie's own widgets, so that it runs in a GLMakie
# window or in a browser (WGLMakie) alike. The caprock in 3D, coloured by
# spill region or depth, with traps and rivers; a click picks a column and
# shows the traps CO2 injected there would fill, in spill order, with a pie
# and a bar chart of their volumes.

# ---------------------------------- Layout ------------------------------------

# Work started by a widget or the mouse runs as a task: in the browser a
# callback holds up the messages the plots and picking need (errors in a task
# are shown, not lost)
_background(f) = @async try
    f()
catch err
    @error "interactive_traps" exception = (err, catch_backtrace())
end

# The figure: a 3D axis over the whole window (a 3D scene is cut off at its
# axis's edge), the title and the panel floating over it. Returns
# (fig, ax, title, panel).
function _viewer_layout(size, view)
    # no padding: the 3D scene reaches the window's edges
    fig = Figure(; size, figure_padding = 0)
    az, el = view
    # :free zooms the whole scene, so close-ups are not clipped at the axis box;
    # shift-click is ours (upstream), not Axis3's rotation reset
    ax = Axis3(fig[1, 1]; zreversed = true, viewmode = :free, clip = false,
               resetrotationkey = false,
               azimuth = deg2rad(az - 90), elevation = deg2rad(el),
               xlabel = "x (m)", ylabel = "y (m)", zlabel = "depth (m)",
               xticks = LinearTicks(3), yticks = LinearTicks(3))
    title = Label(fig[1, 1], ""; font = :bold, valign = :top, padding = (0, 0, 0, 15),
                  tellwidth = false, tellheight = false)
    panel = GridLayout(fig[1, 1]; width = 330, halign = :right, valign = :top,
                       tellwidth = false, tellheight = false, alignmode = Outside(15))
    return (fig, ax, title, panel)
end

# A labelled menu at panel position pos: a plain dropdown (Makie's menus are
# by default also search boxes)
function _menu(pos, label; kw...)
    g = GridLayout(pos)
    Label(g[1, 1], label; halign = :left)
    return Menu(g[2, 1]; searchable = false, kw...)
end

# The panel's widgets below the caller's (`widgets(extra)`, made first: of
# overlapping widgets the first made gets the click, so a caller's dropdown
# list works over ours). Returns a NamedTuple.
function _panel_widgets!(panel, method, exaggeration, toggles, widgets, charts, select)
    # a backdrop under the widgets (made first, so drawn below them)
    backdrop = Box(panel[1:(charts ? 7 : 5), 1]; color = (:white, 0.85),
                   strokevisible = false, alignmode = Outside(-10))
    extra = GridLayout(panel[1, 1])
    # an empty layout has no height, which leaves the panel's undetermined
    Box(extra[1, 1]; height = 0, visible = false)
    own = widgets(extra)
    menus = GridLayout(panel[2, 1])
    m = _menu(menus[1, 1], "Trap analysis"; options = ["node", "cell"],
              default = string(method))
    factors = [1, 2, 5, 10, 20, 50, 100]
    # (an integer default is a position in the options)
    e = _menu(menus[1, 2], "Vertical exaggeration"; options = factors,
              default = something(findfirst(==(exaggeration), factors), 4))

    # on/off switches, two columns of label and toggle (the path's colours only
    # with selection)
    keys = (:traps, :regions, :light, :contours, :red)[1:(select ? 5 : 4)]
    names = ("Traps", "Spill regions", "Lighting", "Contours", "Red path")
    sw = GridLayout(panel[3, 1])
    t = map(enumerate(names[1:length(keys)])) do (k, name)
        r, c = fldmod1(k, 2)
        Label(sw[r, 2c - 1], name; halign = :left)
        Toggle(sw[r, 2c]; active = toggles[k], halign = :left)
    end
    bt = GridLayout(panel[4, 1])
    labels = ("Top view", "Reset view", "Trap detail")[1:(select ? 3 : 2)]
    buttons = [Button(bt[1, k]; label = l) for (k, l) in enumerate(labels)]
    info = Label(panel[5, 1], ""; halign = :left, justification = :left,
                 tellwidth = false)
    # the charts' areas: drawn in pixels (_chart_plots!), as an Axis's own
    # scene does not stay in front of the 3D scene
    area(pos) = Box(pos; height = 140, color = :transparent, strokevisible = false)
    pie, bars = charts ? (area(panel[6, 1]), area(panel[7, 1])) : (nothing, nothing)
    return (backdrop = backdrop, extra = extra, widgets = own, method = m,
            exaggeration = e,
            toggles = NamedTuple{keys}(Tuple(t)), top = buttons[1], reset = buttons[2],
            detail = select ? buttons[3] : nothing, info = info,
            pie = pie, bars = bars)
end

# Blocks drawn in front of the 3D scene (2D blocks are otherwise hidden
# where it is nearer in the depth buffer): near the 2D camera's near limit
# (10 000), the backdrop just behind the widgets on it
function _bring_forward!(blocks, backdrop)
    foreach(b -> translate!(b.blockscene, 0, 0, b === backdrop ? 8500 : 9000), blocks)
end

# Every block in a layout, nested layouts included (not the blocks inside a
# block, e.g. a SliderGrid's: their scenes follow their block's)
_blocks(g::GridLayout) = vcat([_blocks(x) for x in contents(g)]...)
_blocks(b) = [b]

# ----------------------------------- Viewer -----------------------------------

function interactive_traps(s; res = nothing, method::Symbol = :node, traps = true,
                           spill_regions = false, lighting = false, contours = false,
                           path_colours = true, injection_column = 0,
                           exaggeration = 10, view = (-37.5, 30.0), title = "",
                           size = (1700, 1000), widgets = extra -> nothing,
                           charts::Bool = true, selection::Bool = true,
                           on_click = nothing, on_draw = nothing)
    fig, ax, tlabel, panel = _viewer_layout(size, view)
    w = _panel_widgets!(panel, method, exaggeration,
                        (traps, spill_regions, lighting, contours, path_colours),
                        widgets, charts, selection && isnothing(on_click))
    colorbar = _value_colorbar(fig)
    _bring_forward!([tlabel; colorbar; _blocks(panel)], w.backdrop)
    plots = charts ? _chart_plots!(w.pie, w.bars) : (;)
    v = (; fig, ax, title = tlabel, panel, w..., charts = plots, colorbar,
         values = Ref{Any}(nothing), view = view,
         selectable = selection && isnothing(on_click), on_click, on_draw,
         current = Ref{Any}(nothing), selection = Ref((column = 0, upstream = false)),
         plots = Dict{Symbol, Any}(), in_detail = Ref(false), lock = ReentrantLock(),
         pressed_at = Dict{Mouse.Button, Point2f}())
    _connect!(v)
    interactive_traps!(v, s; res, title, injection_column)
    return v
end

function interactive_traps!(v, s; res = nothing, title = "", injection_column = 0)
    v.info.text[] = "Computing the traps ..."
    m = Symbol(v.method.selection[])
    res = something(res, trap_analysis(s; method = m))
    _locked(v) do
        v.current[] = (s = s, res = res, vol = trap_volumes(s, res), title = title)
        # values per column belong to the previous surface
        v.values[] = nothing
        v.colorbar.blockscene.visible[] = false
        v.selection[] = (column = injection_column, upstream = false)
        _draw!(v)
        v.ax.limits[] = (nothing, nothing, nothing)
        reset_limits!(v.ax)
    end
    return v
end

# f() with the viewer's lock held: drawing, recolouring, picking and the trap
# detail run as tasks in a browser, and interleaved they would mix up the
# viewer's plots
_locked(f, v) = lock(f, v.lock)

# The widgets' and the mouse's callbacks
function _connect!(v)
    on(_ -> _background(() -> _reanalyse!(v)), v.method.selection)
    on(_ -> _set_aspect!(v), v.exaggeration.selection)
    for k in (:traps, :regions, :red)
        haskey(v.toggles, k) || continue
        on(_ -> _locked(() -> _update!(v), v), getfield(v.toggles, k).active)
    end
    # lighting is set when a mesh is made: draw anew
    on(_ -> _locked(() -> _draw!(v), v), v.toggles.light.active)
    on(a -> haskey(v.plots, :contours) && (v.plots[:contours].visible[] = a),
       v.toggles.contours.active)
    on(_ -> _reset_view!(v, (0.0, 90.0)), v.top.clicks)
    on(_ -> _reset_view!(v, v.view), v.reset.clicks)
    isnothing(v.detail) || on(_ -> _locked(() -> _toggle_detail!(v), v), v.detail.clicks)
    # the charts follow their areas (moved by a resized window or the text)
    for area in (v.pie, v.bars)
        isnothing(area) && continue
        on(_ -> _locked(() -> _redraw_charts!(v), v), area.layoutobservables.computedbbox)
    end
    on(ev -> _mouse!(v, ev), events(v.ax.scene).mousebutton)
    # over the widgets the axis takes no mouse input (it lies under them all):
    # its turning, panning and zooming switched off there, not during a drag
    e = events(v.fig.scene)
    # a press on the widgets reaches the axis's mouse tracking too, but a menu
    # takes the release: a release sent at the next move ends that "drag"
    pressed_on_widgets = Ref(false)
    on(e.mousebutton; priority = 100) do ev
        ev.action == Mouse.press &&
            (pressed_on_widgets[] = _over_widgets(v, Point2f(e.mouseposition[])))
        return Consume(false)
    end
    on(e.mouseposition; priority = 100) do pos
        isempty(e.mousebuttonstate) || return Consume(false)
        if pressed_on_widgets[]
            pressed_on_widgets[] = false
            empty!(v.pressed_at)
            e.mousebutton[] = Makie.MouseButtonEvent(Mouse.left, Mouse.release)
        end
        _axis_interactions!(v.ax, !_over_widgets(v, Point2f(pos)))
        return Consume(false)
    end
end

# Whether screen position `pos` is on the panel, the colour bar or an open
# dropdown list (which may reach beyond the panel)
function _over_widgets(v, pos)
    pos in v.panel.layoutobservables.computedbbox[] && return true
    v.colorbar.blockscene.visible[] &&
        pos in v.colorbar.layoutobservables.computedbbox[] && return true
    return any(b -> b isa Menu && b.is_open[], _blocks(v.panel))
end

# The axis's mouse interactions (turning, panning, zooming, resets) on or off
function _axis_interactions!(ax, active)
    for (name, (now, _)) in Makie.interactions(ax)
        now == active && continue
        active ? Makie.activate_interaction!(ax, name) :
                 Makie.deactivate_interaction!(ax, name)
    end
end

# The current surface's traps by the method now chosen
function _reanalyse!(v)
    c = v.current[]
    isnothing(c) && return
    column = v.selection[].column
    interactive_traps!(v, c.s; title = c.title, injection_column = column)
end

# The charts drawn anew for the current selection (also in the trap detail,
# whose shorter text moves them)
function _redraw_charts!(v)
    c = v.current[]
    isnothing(c) && return
    _charts!(v, _path(c.res, v.selection[])...)
end

# The view seen from (azimuth, elevation) in co2lab's (MATLAB's) degrees, not
# zoomed or moved
function _reset_view!(v, (az, el))
    v.ax.azimuth[] = deg2rad(az - 90)
    v.ax.elevation[] = deg2rad(el)
    v.ax.zoom_mult[] = 1.0
    v.ax.axis_offset[] = Vec2d(0, 0)
end

# The axis to scale in x and y, its depth exaggerated, for columns `cols`
function _set_aspect!(v, cols = nothing)
    c = v.current[]
    isnothing(c) && return
    s = c.s
    nodes = isnothing(cols) ? eachindex(s.z) : unique(vcat(s.polygons[cols]...))
    span(x) = (r = extrema(view(x, nodes)); max(r[2] - r[1], 1.0))
    m = max(span(s.x), span(s.y))
    v.ax.aspect[] = (span(s.x) / m, span(s.y) / m,
                     v.exaggeration.selection[] * span(s.z) / m)
end

# ---------------------------------- Drawing -----------------------------------

# The caprock with its rivers, outline, faults and depth contours; then the
# colours and the selection (_update!)
function _draw!(v)
    c = v.current[]
    isnothing(c) && return
    s, res = c.s, c.res
    empty!(v.ax)
    empty!(v.plots)
    v.in_detail[] = false
    isnothing(v.detail) || (v.detail.label[] = "Trap detail")
    # one colour per column from the start (a plot of one colour keeps one)
    colours = _column_colours(v, last(_path(res, v.selection[])))
    v.plots[:surface] = _surface_mesh!(v.ax, s; color = colours, offset = 0.0,
                                       colormap = :viridis, smooth = true,
                                       shading = v.toggles.light.active[])
    # the mesh's points' columns, to colour the columns and to pick them
    v.plots[:vcols] = _surface_triangles(s, 0.0)[3]

    δ = _lift(s)
    pts = [Point3f(s.cx[k], s.cy[k], s.zc[k] - δ) for rt in res.rivers for r in rt
           for (a, b) in zip(r[1:end-1], r[2:end]) for k in (a, b)]
    isempty(pts) || linesegments!(v.ax, pts; color = :black, linewidth = 1,
                                  depth_shift = -1f-4)
    plot_ve_outline!(v.ax, s)
    if hasproperty(s, :sealed) && any(s.sealed)
        pts = [Point3f(s.x[n], s.y[n], s.z[n] - δ) for k in findall(s.sealed)
               for n in s.edges[:, k]]
        linesegments!(v.ax, pts; color = :red3, linewidth = 3, depth_shift = -1f-4)
    end
    v.plots[:contours] = _contour_plot!(v.ax, s)
    v.plots[:contours].visible[] = v.toggles.contours.active[]
    _overlay_plots!(v)
    _set_aspect!(v)
    v.title.text[] = c.title
    _update!(v)
    # the caller's own plots, made anew with the caprock's
    isnothing(v.on_draw) || v.on_draw(v)
end

# Ten depth contours on the caprock (co2lab: contourAtlas(data, 10))
function _contour_plot!(ax, s)
    zmin, zmax = extrema(s.z)
    levels = range(zmin, zmax; length = 12)[2:end-1]
    lines = JutulCO2Lab.ve_contour_lines(s, s.z, levels)
    pts = [Point3f(g[i], g[i+1], level - _lift(s)) for (level, segs) in zip(levels, lines)
           for g in segs for i in (1, 3)]
    return linesegments!(ax, pts; color = :black, linewidth = 0.5, depth_shift = -1f-4)
end

# The traps CO2 from the selected column fills, in spill order (trap first,
# then those downstream, or with `upstream` those spilling into it): (trap,
# traps), trap 0 and no traps if the column drains out of the domain
function _path(res, sel)
    sel.column == 0 && return (0, Int[])
    t = res.trap_regions[sel.column]
    t == 0 && return (0, Int[])
    adj = sel.upstream ? permutedims(res.trap_adj) : res.trap_adj
    return (t, downstream_traps((; trap_adj = adj), t))
end

# co2lab's colorizeDepth: the k-th of n traps on a path from bright red to
# nearly black (grey without `red`)
function _path_colour(k, n, red)
    c = max(1 - eps() - (k - 1) / (n + 1), 0.0)
    return red ? RGBAf(0.25 + 0.75c, 0.15, 0.15, 1) : RGBAf(fill(0.25 + 0.5c, 3)..., 1)
end

# A colour per column: depth (viridis) or spill regions, the traps in their
# own colour, the selected path's traps by _path_colour, the selected column
# yellow
function _column_colours(v, path)
    isnothing(v.values[]) || return _value_colours(v)
    c = v.current[]
    s, res = c.s, c.res
    show_traps = v.toggles.traps.active[]
    out = if v.toggles.regions.active[]
        [k > 0 ? _light(_palette_colour(k)) : RGBAf(0.9, 0.9, 0.9, 1)
         for k in res.trap_regions]
    else
        zmin, zmax = extrema(s.zc)
        cmap = Makie.to_colormap(:viridis)
        [RGBAf(Makie.interpolated_getindex(cmap, z, (zmin, max(zmax, zmin + 1))))
         for z in s.zc]
    end
    show_traps && for (k, tr) in enumerate(res.traps)
        tr > 0 && (out[k] = _palette_colour(tr))
    end
    rank = Dict(tr => k for (k, tr) in enumerate(path))
    for (k, tr) in enumerate(res.traps)
        haskey(rank, tr) && (out[k] = _path_colour(rank[tr], length(path),
                                                    v.toggles.red.active[]))
    end
    sel = v.selection[].column
    sel > 0 && (out[sel] = RGBAf(1, 1, 0, 1))
    return out
end

# A colour per column from the values shown (show_column_values!): the
# colormap over their finite range, light grey where NaN, the selected column
# yellow
function _value_colours(v)
    val = v.values[]
    cmap = Makie.to_colormap(val.colormap)
    lims = _value_limits(val.values)
    out = [isfinite(x) ? RGBAf(Makie.interpolated_getindex(cmap, x, lims)) :
           RGBAf(0.85, 0.85, 0.85, 1) for x in val.values]
    sel = v.selection[].column
    sel > 0 && (out[sel] = RGBAf(1, 1, 0, 1))
    return out
end

# The finite range of values (widened if a single value; (0, 1) if none)
function _value_limits(values)
    finite = filter(isfinite, values)
    isempty(finite) && return (0.0, 1.0)
    lo, hi = extrema(finite)
    return hi > lo ? (lo, hi) : (lo - 0.5, hi + 0.5)
end

# The colour bar of values shown, at the bottom left, hidden until used
function _value_colorbar(fig)
    cb = Colorbar(fig[1, 1]; limits = (0.0, 1.0), height = 300, width = 15,
                  halign = :left, valign = :bottom, tellwidth = false,
                  tellheight = false)
    cb.blockscene.visible[] = false
    return cb
end

function show_column_values!(v, values; label = "", colormap = :viridis)
    _locked(v) do
        c = v.current[]
        if isnothing(values)
            v.values[] = nothing
            v.colorbar.blockscene.visible[] = false
        else
            length(values) == length(c.s.zc) ||
                throw(ArgumentError("one value per column needed"))
            v.values[] = (values = collect(Float64, values), label = label,
                          colormap = colormap)
            v.colorbar.limits[] = _value_limits(v.values[].values)
            v.colorbar.colormap[] = colormap
            v.colorbar.label[] = label
            v.colorbar.blockscene.visible[] = true
        end
        _update!(v)
    end
    return v
end

# The colours, the selection's overlays (its traps' outline, rivers, a
# vertical line at the column), the charts and the text
function _update!(v)
    c = v.current[]
    (isnothing(c) || v.in_detail[]) && return
    t, path = _path(c.res, v.selection[])
    colours = _column_colours(v, path)
    v.plots[:surface].color[] = colours[v.plots[:vcols]]
    _selection_overlays!(v, path)
    _charts!(v, t, path)
    v.info.text[] = _report(v, t, path)
end

# The selection's overlays, made empty with each drawing (_draw!) and then
# only given new points: making plots in a browser session is slow
function _overlay_plots!(v)
    lift(colour, width) = linesegments!(v.ax, Point3f[]; color = colour,
                                        linewidth = width, depth_shift = -2f-4)
    v.plots[:outline] = lift(:blue, 2)
    v.plots[:rivers] = lift(:black, 3)
    v.plots[:marker] = lift(:red, 2)
end

# The outline of the path's traps (blue; green upstream), their rivers thick,
# a vertical red line through the selected column
function _selection_overlays!(v, path)
    c = v.current[]
    s, res = c.s, c.res
    δ = _lift(s)
    on_path = in(Set(path))
    inside(k) = k > 0 && on_path(res.traps[k])
    pts = Point3f[Point3f(s.x[n], s.y[n], s.z[n] - 2δ) for k in axes(s.edges, 2)
                  if inside(s.edge_columns[1, k]) != inside(s.edge_columns[2, k])
                  for n in s.edges[:, k]]
    colour = v.selection[].upstream ? :green : :blue
    Makie.update!(v.plots[:outline]; arg1 = pts, color = colour)
    rivers = Point3f[Point3f(s.cx[k], s.cy[k], s.zc[k] - 2δ) for u in path
                     for r in res.rivers[u] for (a, b) in zip(r[1:end-1], r[2:end])
                     for k in (a, b)]
    # (keywords: Makie takes a single positional argument for a Dict)
    Makie.update!(v.plots[:rivers]; arg1 = rivers)
    col = v.selection[].column
    zmin, zmax = extrema(s.z)
    Makie.update!(v.plots[:marker]; arg1 = col == 0 ? Point3f[] :
                  [Point3f(s.cx[col], s.cy[col], z) for z in (zmin, zmax)])
end

# co2lab's pie (the trap, the traps further on, the other traps' volume) and
# its bars (each trap's volume, in spill order, on a log scale), drawn in
# pixels in the areas `pie` and `bars`: made once per viewer, then only given
# new points (making plots in a browser is slow). Returns a NamedTuple of the
# plots.
function _chart_plots!(pie, bars)
    sp, sb = pie.blockscene, bars.blockscene
    blank = [Point2f[(0, 0), (0, 0), (0, 0)] for _ in 1:3]
    wedges = poly!(sp, blank; color = _PIE_COLOURS, strokecolor = :white,
                   strokewidth = 1)
    title = text!(sp, Point2f(0); text = "Volumes (m³)", font = :bold, fontsize = 14,
                  align = (:center, :top))
    names = text!(sp, fill(Point2f(0), 3); text = fill("", 3), color = _PIE_COLOURS,
                  fontsize = 12, align = (:left, :center))
    rects = poly!(sb, [Rect2f(0, 0, 1, 1)]; color = [RGBAf(1, 0, 0, 1)])
    frame = lines!(sb, Point2f[]; color = :black, linewidth = 1)
    ticks = text!(sb, fill(Point2f(0), 3); text = fill("", 3), fontsize = 11,
                  align = (:right, :center))
    labels = text!(sb, fill(Point2f(0), 2); text = ["Trapped volume (m³)",
                                                    "Traps in spill order"],
                   fontsize = 12, rotation = [pi / 2, 0.0],
                   align = [(:center, :bottom), (:center, :top)])
    return (wedges = wedges, title = title, names = names, rects = rects, frame = frame,
            ticks = ticks, labels = labels)
end

# The pie's colours: the trap, the traps further on, the other traps
const _PIE_COLOURS = [RGBf(0, 0, 0.7), RGBf(0.1, 0.8, 0.1), RGBf(0.6, 0, 0)]

# The charts for trap t and its path (hidden for none)
function _charts!(v, t, path)
    isempty(v.charts) && return
    foreach(p -> p.visible[] = t > 0, values(v.charts))
    t == 0 && return
    vol = v.current[].vol
    primary, further = vol[t], sum(vol[path[2:end]]; init = 0.0)
    # (co2lab: zero shown as 1 m³, for the pie)
    rest = max(sum(vol) - primary - further, eps())
    other = v.selection[].upstream ? "Upstream" : "Migration"
    _pie!(v, [primary, max(further, 1.0), rest], ["Primary", other, "Not filled"])
    n = length(path)
    _bars!(v, max.(vol[path], 1.0), [_path_colour(k, n, v.toggles.red.active[])
                                     for k in 1:n])
end

# The pie in its area: wedges of `values`, the title above, `names` beside
function _pie!(v, values, names)
    r = v.pie.layoutobservables.computedbbox[]
    (x0, y0), (w, h) = r.origin, r.widths
    R = min(0.3w, 0.5h - 14)
    c = Point2f(x0 + 0.35w, y0 + 0.5h - 8)
    ends = 2pi .* cumsum(values) ./ sum(values)
    arc(a, b) = [c .+ R .* Point2f(cos(x), sin(x)) for x in range(a, b; length = 40)]
    wedge(a, b) = pushfirst!(arc(a, b), c)
    Makie.update!(v.charts.wedges; arg1 = [wedge(a, b) for (a, b) in
                                           zip([0.0; ends[1:end-1]], ends)])
    # the title at the top, the names beside the pie
    Makie.update!(v.charts.title; arg1 = Point2f(x0 + 0.5w, y0 + h))
    at = [Point2f(x0 + 0.72w, c[2] + 18 - 18k) for k in 0:2]
    Makie.update!(v.charts.names; arg1 = at, text = names)
end

# The bars in their area: `values` on a log scale from 1 m³, with a frame,
# decade marks and the axis labels
function _bars!(v, values, colours)
    r = v.bars.layoutobservables.computedbbox[]
    (x0, y0), (w, h) = r.origin, r.widths
    # the plot box inside the area, room left for the marks and labels
    bx, by, bw, bh = x0 + 62, y0 + 20, w - 66, h - 26
    top = max(1, ceil(log10(maximum(values))))
    n = length(values)
    dx = bw / n
    gap = n > 150 ? 0.0 : 0.2dx
    rects = [Rect2f(bx + (k - 1) * dx + gap / 2, by, dx - gap, bh * log10(x) / top)
             for (k, x) in enumerate(values)]
    Makie.update!(v.charts.rects; arg1 = rects, color = colours)
    corners = Point2f[(bx, by), (bx + bw, by), (bx + bw, by + bh), (bx, by + bh),
                      (bx, by)]
    Makie.update!(v.charts.frame; arg1 = corners)
    decades = round.(Int, range(0, top; length = 3))
    marks = [Point2f(bx - 4, by + bh * d / top) for d in decades]
    Makie.update!(v.charts.ticks; arg1 = marks,
                  text = [d == 0 ? "1" : "1e$d" for d in decades])
    Makie.update!(v.charts.labels; arg1 = [Point2f(x0 + 12, by + bh / 2),
                                           Point2f(bx + bw / 2, by - 4)])
end

# The text under the switches: the selected column and its traps
function _report(v, t, path)
    c = v.current[]
    s, vol = c.s, c.vol
    head = @sprintf("%d columns, %d traps, %.3g Mm³ in all", length(s.zc), length(vol),
                    sum(vol; init = 0.0) / 1e6)
    col = v.selection[].column
    v.selectable || return head
    col == 0 && return head * "\n\nClick a column (shift: upstream)."
    at = @sprintf("\n\nColumn %d at (%.0f, %.0f), %.0f m", col, s.cx[col], s.cy[col],
                     s.zc[col])
    t == 0 && return head * at * "\nspills out of the domain"
    dir = v.selection[].upstream ? "spill into it" : "it spills into"
    trap = @sprintf("\nTrap %d: %.3g Mm³, spill depth %.0f m", t, vol[t] / 1e6,
                    c.res.trap_z[t])
    further = @sprintf("\n%d trap(s) %s: %.3g Mm³", length(path) - 1, dir,
                       sum(vol[path[2:end]]; init = 0.0) / 1e6)
    return head * at * trap * further
end

# ----------------------------------- Mouse ------------------------------------

# A click is a press and release without moving (a drag turns or pans the
# view): left picks a column (shift-left or middle: upstream, as co2lab),
# right clears the selection; clicks on the panel are the widgets'
function _mouse!(v, ev)
    (v.selectable || !isnothing(v.on_click)) || return
    ev.button in (Mouse.left, Mouse.middle, Mouse.right) || return
    scene = v.ax.scene
    pos = Point2f(events(scene).mouseposition[])
    if ev.action == Mouse.press
        v.pressed_at[ev.button] = pos
        return
    end
    sum(abs2, pos .- get(v.pressed_at, ev.button, Point2f(Inf))) < 9 || return
    (isnothing(v.current[]) || v.in_detail[]) && return
    _over_widgets(v, pos) && return
    # a caller's clicks: the column (0: none) and the button, to on_click
    if !isnothing(v.on_click)
        origin = minimum(v.ax.scene.viewport[])
        _background(() -> _locked(() -> v.on_click(_pick_column(v, pos .- origin),
                                                     ev.button), v))
        return
    end
    if ev.button == Mouse.right
        _locked(v) do
            v.selection[] = (column = 0, upstream = false)
            _update!(v)
        end
        return
    end
    keys = events(scene).keyboardstate
    upstream = ev.button == Mouse.middle || Keyboard.left_shift in keys ||
               Keyboard.right_shift in keys
    _background(() -> _select!(v, pos .- minimum(v.ax.scene.viewport[]), upstream))
end

# The column under screen position `pos` (pixels in the axis's scene), and
# selected
function _select!(v, pos, upstream)
    _locked(v) do
        col = _pick_column(v, pos)
        col == 0 && return
        v.selection[] = (column = col, upstream = upstream)
        _update!(v)
    end
end

# The caprock column drawn at screen position `pos` (0: none): the nearest of
# the surface's triangles that contain it on screen, projected with the
# camera's matrices. Done here, not by the backend's picking: WGLMakie's
# picks a stale view after some camera changes (e.g. a view from above), and
# lines drawn over the caprock do not get in the way.
function _pick_column(v, pos)
    plot = v.plots[:surface]
    m = plot[1][]
    M = v.ax.scene.camera.projectionview[] * plot.model[]
    w = widths(v.ax.scene.viewport[])
    # each point on screen (pixels) and its depth (smaller: nearer); NaN behind
    screen = map(Makie.GeometryBasics.coordinates(m)) do p
        c = M * Vec4f(p[1], p[2], p[3], 1)
        c[4] <= 0 && return Point3f(NaN)
        Point3f((c[1] / c[4] + 1) * w[1] / 2, (c[2] / c[4] + 1) * w[2] / 2, c[3] / c[4])
    end
    best, depth = 0, Inf
    for f in Makie.GeometryBasics.faces(m)
        a, b, c = (screen[convert(Int, k)] for k in f)
        d = _depth_in_triangle(pos, a, b, c)
        d < depth && ((best, depth) = (convert(Int, f[1]), d))
    end
    return best == 0 ? 0 : v.plots[:vcols][best]
end

# The depth at screen point q in the screen triangle (a, b, c) (Inf outside)
function _depth_in_triangle(q, a, b, c)
    det = (b[2] - c[2]) * (a[1] - c[1]) + (c[1] - b[1]) * (a[2] - c[2])
    (isnan(det) || det == 0) && return Inf
    la = ((b[2] - c[2]) * (q[1] - c[1]) + (c[1] - b[1]) * (q[2] - c[2])) / det
    lb = ((c[2] - a[2]) * (q[1] - c[1]) + (a[1] - c[1]) * (q[2] - c[2])) / det
    lc = 1 - la - lb
    min(la, lb, lc) < 0 && return Inf
    return la * a[3] + lb * b[3] + lc * c[3]
end

# -------------------------------- Trap detail ---------------------------------

# co2lab's trap figure (its right click), here in the same axis: the
# selected trap's floor (filled to the spill depth) red, its ceiling blue, the
# caprock around it faint; the button again returns. The caprock's plots are
# only hidden meanwhile: made anew, they would be placed by the trap's axis
# limits in a browser (and picking would miss), and making them is slow there.
function _toggle_detail!(v)
    c = v.current[]
    isnothing(c) && return
    v.in_detail[] && return _leave_detail!(v)
    t, _ = _path(c.res, v.selection[])
    t == 0 && (v.info.text[] = "Select a trap first."; return)
    v.in_detail[] = true
    v.detail.label[] = "Back"
    foreach(p -> p.visible[] = false, v.ax.scene.plots)
    n = length(v.ax.scene.plots)
    _detail_plots!(v, t)
    v.plots[:detail] = v.ax.scene.plots[n+1:end]
end

# The trap's plots removed, the caprock's shown again, the axis fitted to it
function _leave_detail!(v)
    foreach(p -> delete!(v.ax, p), pop!(v.plots, :detail))
    foreach(p -> p.visible[] = true, v.ax.scene.plots)
    v.plots[:contours].visible[] = v.toggles.contours.active[]
    v.in_detail[] = false
    v.detail.label[] = "Trap detail"
    v.ax.limits[] = (nothing, nothing, nothing)
    reset_limits!(v.ax)
    _set_aspect!(v)
    _update!(v)
end

# The trap's plots, the axis fitted to its surroundings
function _detail_plots!(v, t)
    s, res = v.current[].s, v.current[].res
    trap = findall(==(t), res.traps)
    lo = [minimum(s.cx[trap]), minimum(s.cy[trap])]
    hi = [maximum(s.cx[trap]), maximum(s.cy[trap])]
    # 30% around it, and at least two columns (a trap may be one column)
    r = max.(0.3 .* (hi .- lo), 2sqrt(maximum(s.area[trap])))
    near = findall(k -> all(lo .- r .<= (s.cx[k], s.cy[k]) .<= hi .+ r), eachindex(s.zc))
    around = setdiff(near, trap)
    sub(s, cols) = (; s..., polygons = s.polygons[cols])

    plot_ve_surface!(v.ax, sub(flatten_traps(s, res), trap); color = :red,
                     strokewidth = 0, shading = true)
    plot_ve_surface!(v.ax, sub(s, trap); color = (:blue, 0.8), strokewidth = 0,
                     shading = true)
    # the caprock around, coloured by depth, faint
    cmap, zr = Makie.to_colormap(:viridis), extrema(s.zc[near])
    faint = map(around) do k
        c = Makie.interpolated_getindex(cmap, s.zc[k], zr)
        RGBAf(c.r, c.g, c.b, 0.3)
    end
    plot_ve_surface!(v.ax, sub(s, around); color = faint, strokewidth = 0, shading = true)
    nodes = unique(vcat(s.polygons[near]...))
    ext(x) = extrema(view(x, nodes))
    v.ax.limits[] = (ext(s.x)..., ext(s.y)..., ext(s.z)...)
    _set_aspect!(v, near)
    v.info.text[] = "Trap $t: floor (filled to its spill depth) red,\n" *
                    @sprintf("ceiling blue; %.3g Mm³", v.current[].vol[t] / 1e6)
end
