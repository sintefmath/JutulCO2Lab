# The trapping inventory over a run (co2lab's plotTrappingDistribution).

# co2lab's getInventoryColors, as plotTrappingDistribution assigns them
function _inventory_colours()
    mix(a, b, c) = RGBf((c .* a .+ (100 - c) .* b) ./ 100...)
    red, green, yellow = (1, 0, 0), (0, 1, 0), (1, 1, 0)
    orange, white, black = (1, 0.65, 0), (1, 1, 1), (0, 0, 0)
    c = [mix(green, black, 50), mix(green, black, 75), green, mix(green, white, 50),
         mix(yellow, white, 80), mix(orange, white, 80), mix(orange, red, 30)]
    # subscale: halfway between residual in plume and structural plume
    subscale = RGBf(((c[5].r + c[4].r) / 2, (c[5].g + c[4].g) / 2,
                     (c[5].b + c[4].b) / 2)...)
    # dissolved, structural residual, residual, residual in plume, subscale,
    # structural plume, free plume, exited
    return [c[1], c[2], c[3], c[4], subscale, c[5], c[6], c[7]]
end

const _INVENTORY_LABELS = ["Dissolved", "Structural residual", "Residual",
                           "Residual in plume", "Structural subscale",
                           "Structural plume", "Free plume", "Exited"]

function plot_trapping_history(history; title = "", kwarg...)
    fig = Figure(size = (1000, 550))
    ax = Axis(fig[1, 1]; title = title, xlabel = "years since simulation start",
              ylabel = "mass (Mt)")
    plot_trapping_history!(ax, history; kwarg...)
    return fig
end

function plot_trapping_history!(ax, history; legend = true)
    years = history.t ./ (365.25 * 86400)
    # the exited mass on top (negative round-off clipped), in Mt
    masses = [[collect(values(inv)); max(e, 0.0)] ./ 1e9
              for (inv, e) in zip(history.inventory, history.exited)]
    colours = _inventory_colours()
    # stacked bands: each category from the top of the previous one
    lower = zeros(length(years))
    for k in eachindex(_INVENTORY_LABELS)
        upper = lower .+ getindex.(masses, k)
        # categories that stay zero get no legend entry
        label = any(>(0), getindex.(masses, k)) ? _INVENTORY_LABELS[k] : nothing
        band!(ax, years, lower, upper; color = colours[k], label = label)
        lower = upper
    end
    legend && axislegend(ax; position = :lt)
    return ax
end
