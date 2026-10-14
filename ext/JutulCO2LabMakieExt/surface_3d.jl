# The top surface of a VE domain in 3D, at the caprock's depths, as co2lab draws
# its top-surface grid (plotCellData, plotGrid on Gt), and the Axis3 methods of
# the map, outline and trap plots. Public functions are documented in
# src/plotting.jl.

# ---------------------------------- The axis ----------------------------------

# Depth down, x and y to scale, MATLAB's view angles
function ve_axis3(pos, surface; exaggeration = 10.0, view = (-37.5, 30.0), title = "")
    dx, dy, dz = (-(reverse(extrema(v))...) for v in (surface.x, surface.y, surface.z))
    # x and y to scale, depth exaggerated
    m = max(dx, dy)
    # MATLAB's view(az, el) looks from azimuth az measured from −y
    az, el = view
    # seen from straight above, the depth axis only clutters the view
    top = el >= 89
    return Axis3(pos; title = title, zreversed = true,
                 aspect = (dx / m, dy / m, exaggeration * max(dz, eps()) / m),
                 azimuth = deg2rad(az - 90), elevation = deg2rad(el),
                 xlabel = "x (m)", ylabel = "y (m)", zlabel = "depth (m)",
                 xticks = LinearTicks(3), yticks = LinearTicks(3),
                 zticks = LinearTicks(3), zticklabelsvisible = !top,
                 zlabelvisible = !top, zticksvisible = !top)
end

# ------------------------------- Surface in 3D --------------------------------

# The columns' top faces as one mesh, and their edges
function plot_ve_surface!(ax, surface; color = (:steelblue, 0.6), offset = 0.0,
                          colormap = :viridis, shading = false, strokecolor = :black,
                          strokewidth = 1, smooth = true)
    _surface_mesh!(ax, surface; color, offset, colormap, shading, smooth)
    if strokewidth > 0
        x, y, z = surface.x, surface.y, surface.z .+ offset
        edges = Point3f[Point3f(x[n], y[n], z[n]) for p in surface.polygons
                 for (k, a) in enumerate(p) for n in (a, p[mod1(k + 1, length(p))])]
        linesegments!(ax, edges; color = strokecolor, linewidth = strokewidth,
                      transparency = _translucent(strokecolor))
    end
    return ax
end

# The surface's columns as one mesh (plot_ve_surface! without the edges).
# Returns the mesh plot. Modifies: `ax`.
function _surface_mesh!(ax, surface; color, offset, colormap, shading, smooth)
    pts, tris, col, node = _surface_triangles(surface, offset)
    # one value or colour per column, or one colour for all
    c = color isa AbstractVector ? color[col] : color
    # smooth: the light blends across column edges (normals shared at the
    # nodes); otherwise Makie lights each column as a facet of its own
    m = smooth ? Makie.GeometryBasics.Mesh(pts, tris;
                     normal = _node_normals(surface, pts, tris, node)) :
                 Makie.GeometryBasics.Mesh(pts, tris)
    # translucent colours must be flagged: GLMakie otherwise writes them to its
    # depth buffer, hiding what lies behind them
    return mesh!(ax, m; color = c, colormap = colormap, shading = shading,
                 transparency = _translucent(color))
end

# Whether a colour (or any colour of a vector of colours) is translucent;
# values drawn through a colormap are not
_translucent(c::AbstractVector{<:Real}) = false
_translucent(c::AbstractVector) = any(_translucent, c)
_translucent(c) = Makie.to_color(c).alpha < 1

# Each polygon as a fan of triangles around its node mean: points, triangles,
# the column of each point and its node (0: the polygon's centre). The
# triangles turn the same way in map view whatever the grid's (so that the lit
# side is the upper one).
function _surface_triangles(surface, offset)
    x, y, z = surface.x, surface.y, surface.z .+ offset
    Face = Makie.GeometryBasics.GLTriangleFace
    pts, tris, col, node = Point3f[], Face[], Int[], Int[]
    for (c, p) in enumerate(surface.polygons)
        q = unique(p)
        length(q) < 3 && continue
        k0 = length(pts) + 1
        push!(pts, Point3f(sum(x[q]) / length(q), sum(y[q]) / length(q),
                           sum(z[q]) / length(q)))
        append!(pts, Point3f.(x[q], y[q], z[q]))
        append!(col, fill(c, length(q) + 1))
        append!(node, [0; q])
        # clockwise in map view: with depth as z, the normal points up
        ccw = _signed_area(x[q], y[q]) > 0
        for k in eachindex(q)
            a, b = k0 + k, k0 + mod1(k + 1, length(q))
            push!(tris, ccw ? Face(k0, b, a) : Face(k0, a, b))
        end
    end
    return (pts, tris, col, node)
end

# Each triangle's (unnormalised) normal, from the surface's Float64 positions
# (Float32 map coordinates of ~1e6 m lose the small depth differences)
function _triangle_normals(surface, pts, tris, node)
    p(k) = node[k] == 0 ? Vec3{Float64}(pts[k]...) :
           Vec3(surface.x[node[k]], surface.y[node[k]], surface.z[node[k]])
    normal(a, b, c) = cross(p(b) - p(a), p(c) - p(a))
    return [normal(convert.(Int, Tuple(f))...) for f in tris]
end

# Normals per point: at a node the mean over every triangle touching the node,
# in any column; at a polygon's centre over its own triangles
function _node_normals(surface, pts, tris, node)
    tn = _triangle_normals(surface, pts, tris, node)
    at_node = fill(Vec3(0.0, 0.0, 0.0), length(surface.z))
    at_point = fill(Vec3(0.0, 0.0, 0.0), length(pts))
    for (f, n) in zip(tris, tn), k in convert.(Int, Tuple(f))
        node[k] == 0 ? (at_point[k] += n) : (at_node[node[k]] += n)
    end
    nrm(v) = Vec3f(v / max(norm(v), eps()))
    return [nrm(node[k] == 0 ? at_point[k] : at_node[node[k]]) for k in eachindex(pts)]
end

# Twice the signed area of a polygon in map view (positive: anticlockwise)
function _signed_area(x, y)
    n = length(x)
    return sum(x[k] * y[mod1(k + 1, n)] - x[mod1(k + 1, n)] * y[k] for k in 1:n)
end

# A small lift (towards shallower) to draw lines and overlays over the surface
_lift(surface) = 2e-3 * max(-(reverse(extrema(surface.z))...), 1.0)

# The lateral boundary just above the caprock
function plot_ve_outline!(ax::Axis3, surface; color = :black, linewidth = 2,
                          label = nothing)
    x, y, z = surface.x, surface.y, surface.z .- _lift(surface)
    pts = Point3f[Point3f(x[n], y[n], z[n]) for e in surface.outline for n in e]
    linesegments!(ax, pts; color = color, linewidth = linewidth, label = label,
                  depth_shift = -1f-4)
    return ax
end

# ------------------------------------ Maps ------------------------------------

# The map on the caprock (no plumes)
function plot_ve_map!(ax::Axis3, surface; background = nothing, colormap = :viridis,
                      depth_contours = 0, plumes = [], wells = Int[], regions = nothing,
                      region_color = (:red, 0.3), outline = true, kwarg...)
    isempty(plumes) || throw(ArgumentError("plumes are drawn on map views (Axis) only"))
    δ = _lift(surface)
    isnothing(background) ||
        plot_ve_surface!(ax, surface; color = collect(Float64, background),
                         colormap = colormap, shading = true, strokewidth = 0)

    # depth contours at their depth
    if depth_contours > 0
        zmin, zmax = extrema(surface.z)
        levels = range(zmin, zmax; length = depth_contours + 2)[2:end-1]
        lines = JutulCO2Lab.ve_contour_lines(surface, surface.z, levels)
        pts = Point3f[Point3f(s[i], s[i+1], level - δ)
                      for (level, segs) in zip(levels, lines)
                      for s in segs for i in (1, 3)]
        linesegments!(ax, pts; color = :black, linewidth = 0.5, depth_shift = -1f-4)
    end

    # regions (e.g. traps) as an overlay just above the surface
    if !isnothing(regions)
        part = (; surface..., polygons = surface.polygons[findall(Bool.(regions))])
        plot_ve_surface!(ax, part; color = region_color, offset = -δ, strokewidth = 0)
    end

    # wells as vertical lines above their columns (from the polygon's mean depth)
    if !isempty(wells)
        ztop(c) = (q = unique(surface.polygons[c]); sum(surface.z[q]) / length(q))
        pts = Point3f[Point3f(surface.xc[c], surface.yc[c], ztop(c) - t)
               for c in wells for t in (0.0, 50δ)]
        linesegments!(ax, pts; color = :black, linewidth = 3)
    end
    outline && plot_ve_outline!(ax, surface)
    return ax
end

# ------------------------------------ Traps -----------------------------------

# The trap colours on the lit caprock, with rivers, labels, outline, faults
function plot_traps!(ax::Axis3, s, res; regions = true, trap_values = nothing,
                     colormap = :viridis, rivers = true, labels = false,
                     outline = true, faults = true)
    colours = _trap_colours(res; regions, trap_values, colormap)
    plot_ve_surface!(ax, s; color = colours, shading = true, strokewidth = 0)
    δ = _lift(s)

    # rivers through the column centroids, at the caprock
    if rivers
        # typed: a surface whose traps have no rivers gives an empty list
        pts = Point3f[Point3f(s.cx[c], s.cy[c], s.zc[c] - δ) for rt in res.rivers
                      for r in rt for (a, b) in zip(r[1:end-1], r[2:end])
                      for c in (a, b)]
        linesegments!(ax, pts; color = :black, linewidth = 1, depth_shift = -1f-4)
    end
    if labels
        # one label per trap, at its first top column
        top = unique(c -> res.traps[c], res.top)
        text!(ax, [Point3f(s.cx[c], s.cy[c], s.zc[c] - 5δ) for c in top];
              text = string.(res.traps[top]), align = (:center, :center),
              fontsize = 11)
    end
    outline && plot_ve_outline!(ax, s)
    if faults && hasproperty(s, :sealed) && any(s.sealed)
        pts = Point3f[Point3f(s.x[n], s.y[n], s.z[n] - δ) for k in findall(s.sealed)
               for n in s.edges[:, k]]
        linesegments!(ax, pts; color = :red3, linewidth = 3, depth_shift = -1f-4)
    end
    return ax
end

# Trap k's colour (tab20, repeating), and a light version for its spill region
const _TAB20 = Makie.to_colormap(:tab20)
_palette_colour(k) = RGBAf(_TAB20[mod1(k, length(_TAB20))])
_light(c) = RGBAf(0.45c.r + 0.55, 0.45c.g + 0.55, 0.45c.b + 0.55, 1)

# A colour per column: spill regions light (grey where draining out; neutral
# without `regions`), traps in their region's full colour or by `trap_values`
function _trap_colours(res; regions, trap_values, colormap)
    colour, light = _palette_colour, _light
    grey = RGBAf(0.9, 0.9, 0.9, 1)
    out = [regions && k > 0 ? light(colour(k)) : grey for k in res.trap_regions]

    cmap = Makie.to_colormap(colormap)
    limits = isnothing(trap_values) ? (0.0, 1.0) : Float64.(extrema(trap_values))
    for (c, t) in enumerate(res.traps)
        t == 0 && continue
        out[c] = isnothing(trap_values) ? colour(t) :
                 RGBAf(Makie.interpolated_getindex(cmap, trap_values[t], limits))
    end
    return out
end
