# Map view of a VE domain (co2lab's mapPlot) on the columns' top faces.

# ---------------------------------- The map -----------------------------------

function plot_ve_map(surface; title = "", colorbar = true, kwarg...)
    fig, ax = _map_figure(title)
    plot_ve_map!(ax, surface; kwarg...)
    p = _background_plot(ax)
    colorbar && !isnothing(p) && Colorbar(fig[1, 2], p)
    return fig
end

function plot_ve_map!(ax, surface; background = nothing, colormap = :viridis,
                      depth_contours = 40, plumes = [], plume_threshold = 0.3,
                      plume_colors = [:red, :blue, :cyan, :green, :magenta],
                      labels = nothing, wells = Int[], regions = nothing,
                      region_color = (:red, 0.3), outline = true)
    polys = _polygons(surface)
    if !isnothing(background)
        poly!(ax, polys; color = collect(Float64, background), colormap = colormap,
              strokewidth = 0)
    end

    # caprock depth contours on the surface's node depths
    if depth_contours > 0
        zmin, zmax = extrema(surface.z)
        levels = range(zmin, zmax; length = depth_contours + 2)[2:end-1]
        segments = reduce(vcat, JutulCO2Lab.ve_contour_lines(surface, surface.z, levels))
        linesegments!(ax, _points(segments); color = :black, linewidth = 0.5)
    end

    # regions (e.g. traps) as translucent cells
    if !isnothing(regions)
        poly!(ax, polys[findall(Bool.(regions))]; color = region_color, strokewidth = 0)
    end

    # plume outlines where each field crosses the threshold
    length(plumes) <= length(plume_colors) ||
        throw(ArgumentError("more plumes than plume_colors"))
    for (k, field) in enumerate(plumes)
        s = only(JutulCO2Lab.ve_contour_lines(surface, field, [plume_threshold]))
        label = isnothing(labels) ? nothing : labels[k]
        linesegments!(ax, _points(s); color = plume_colors[k], linewidth = 2,
                      label = label)
    end
    isnothing(labels) || isempty(plumes) || axislegend(ax; position = :rb)

    if !isempty(wells)
        scatter!(ax, surface.xc[wells], surface.yc[wells]; color = :black,
                 strokecolor = :yellow, strokewidth = 1, markersize = 12)
    end
    outline && plot_ve_outline!(ax, surface)
    return ax
end

function plot_ve_outline!(ax, surface; color = :black, linewidth = 2, label = nothing)
    x, y = surface.x, surface.y
    pts = [Point2f(x[n], y[n]) for e in surface.outline for n in e]
    linesegments!(ax, pts; color = color, linewidth = linewidth, label = label)
    return ax
end

# ---------------------------------- Helpers -----------------------------------

# Column polygons as Makie points
_polygons(s) = [Point2f.(s.x[p], s.y[p]) for p in s.polygons]

# Segments (x1, y1, x2, y2) as point pairs for linesegments!
_points(segments) = [Point2f(p...) for s in segments
                     for p in ((s[1], s[2]), (s[3], s[4]))]

# The background poly plot of a map axis, if any (for its colour bar)
function _background_plot(ax)
    is_field(p) = p isa Poly && p.color[] isa AbstractVector{<:Real}
    i = findfirst(is_field, ax.scene.plots)
    return isnothing(i) ? nothing : ax.scene.plots[i]
end

# ------------------------------- Surface in 3D --------------------------------

function plot_ve_surface!(ax, surface; color = (:steelblue, 0.6), offset = 0.0,
                          strokecolor = :black, strokewidth = 1)
    x, y, z = surface.x, surface.y, surface.z .+ offset
    Face = Makie.GeometryBasics.GLTriangleFace
    pts, tris = Point3f[], Face[]
    # each polygon as a fan around its node mean
    for p in surface.polygons
        q = unique(p)
        length(q) < 3 && continue
        c = length(pts) + 1
        push!(pts, Point3f(sum(x[q]) / length(q), sum(y[q]) / length(q),
                           sum(z[q]) / length(q)))
        append!(pts, Point3f.(x[q], y[q], z[q]))
        for k in eachindex(q)
            push!(tris, Face(c, c + k, c + (k == length(q) ? 1 : k + 1)))
        end
    end
    mesh!(ax, pts, tris; color = color, shading = NoShading)

    edges = [Point3f(x[n], y[n], z[n]) for p in surface.polygons
             for (k, a) in enumerate(p) for n in (a, p[k == length(p) ? 1 : k + 1])]
    linesegments!(ax, edges; color = strokecolor, linewidth = strokewidth)
    return ax
end
