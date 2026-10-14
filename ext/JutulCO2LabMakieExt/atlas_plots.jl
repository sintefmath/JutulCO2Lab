# CO2 Storage Atlas data in map view: raw maps (co2lab's contourAtlas), the
# coastline and the shelf's wells. Public functions are documented in
# src/plotting.jl.

# A new figure with one map axis
function plot_atlas_map(map; title = "$(map.name) $(map.variant)", kwarg...)
    fig, ax = _map_figure(title)
    plot_atlas_map!(ax, map; kwarg...)
    # a heatmap just drawn gets a colour bar; contour lines none
    p = ax.scene.plots[end]
    p isa Heatmap && Colorbar(fig[1, 2], p)
    return fig
end

# A heatmap or contour lines (2D axis)
function plot_atlas_map!(ax, map; style = :heatmap, levels = 10, colormap = :viridis)
    x, y = JutulCO2Lab._map_positions(map.meta)
    if style == :heatmap
        heatmap!(ax, x, y, map.data; colormap = colormap, nan_color = :transparent)
    elseif style == :contour
        # Makie's contour fails on NaN, which fills every map outside its
        # formation: trace the lines on the lattice instead
        surface = JutulCO2Lab._atlas_map_surface(map)
        vmin, vmax = extrema(filter(!isnan, map.data))
        lv = levels isa Integer ?
             range(vmin, vmax; length = levels + 2)[2:end-1] : levels
        lines = JutulCO2Lab.ve_contour_lines(surface, surface.z, lv)
        colors = cgrad(colormap)
        for (l, s) in zip(lv, lines)
            c = colors[(l - first(lv)) / max(last(lv) - first(lv), eps())]
            linesegments!(ax, _points(s); color = c, linewidth = 1)
        end
    else
        throw(ArgumentError("style must be :heatmap or :contour, got :$style"))
    end
    return ax
end

# In 3D: the map as a surface with its values as height (co2lab's surf)
function plot_atlas_map!(ax::Axis3, map; colormap = :viridis, kwarg...)
    x, y = JutulCO2Lab._map_positions(map.meta)
    surface!(ax, x, y, map.data; colormap = colormap, nan_color = :transparent)
    return ax
end

# Coastline and wells for orientation
function plot_atlas_context!(ax; coastline = true, wells = true)
    if coastline
        for l in atlas_coastline()
            lines!(ax, l[:, 1], l[:, 2]; color = :gray40, linewidth = 1)
        end
    end
    if wells
        w = atlas_wells()
        scatter!(ax, w.x, w.y; color = :black, markersize = 2)
    end
    return ax
end
