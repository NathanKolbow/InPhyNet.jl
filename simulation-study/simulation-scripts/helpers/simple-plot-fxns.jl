using Plots, StatsPlots


function plot_negative_one_prop(netid, std0)
    stepsize = 0.1
    p = plot(title="Network: $(netid)", xlabel="std0 multiplier", ylabel="Proportion of runs that succeeded")

    startcenter = minimum(gausserrors) + stepsize / 2
    endcenter = maximum(gausserrors) - stepsize / 2
    first_non_one = -1.
    start_all_zero = -1.
    last_seen_zero = true

    for x = startcenter:stepsize:endcenter
        iterruns = (gausserrors .>= x - stepsize / 2) .&& (gausserrors .< x + stepsize / 2)
        if sum(iterruns) == 0 continue end
        yval = mean(esterrors[iterruns] .!= -1.)
        scatter!([x / std0], [yval], label=nothing, color="blue")

        if yval < 0.99 && first_non_one == -1.
            first_non_one = x / std0
        end
        if yval < 0.01 && !last_seen_zero
            start_all_zero = x / std0
            last_seen_zero = true
        elseif yval > 0.01 && last_seen_zero
            last_seen_zero = false
        end
    end
    plot!([first_non_one - stepsize / 8, first_non_one - stepsize / 8], [0, 1], color="red", label=nothing)
    plot!([start_all_zero - stepsize / 8, start_all_zero - stepsize / 8], [0, 1], color="red", label=nothing)
    display(p)

    # Print some relevant values
    
    println("$(netid), ($(round(first_non_one, digits=2)), $(round(start_all_zero, digits=2))) - prop-1: $(round(mean(esterrors .== -1), digits=2))")

    return p
end
plot_negative_one_prop() = plot_negative_one_prop(netid, upperTriangStd(D))


"""
x-axis: grouped gaussian errors
y-axis: sum of constraint diffs boxplots
"""
function plot_negative_one_dists(netid, std0)
    neg_one = esterrors .== -1.

    stepsize = 0.1
    p = plot(title="Network: $(netid)", xlabel="std0 multiplier", ylabel="Sum of constraint diffs")

    startcenter = minimum(gausserrors) + stepsize / 2
    endcenter = maximum(gausserrors) - stepsize / 2

    for x = startcenter:stepsize:endcenter
        iterruns = (gausserrors .>= x - stepsize / 2) .&& (gausserrors .< x + stepsize / 2)
        if sum(neg_one .&& iterruns) == 0
            scatter!([x], [0], label=nothing, color="blue")
        else
            boxplot!([x / std0], constraintdiffs[neg_one .&& iterruns], label=nothing, color="red", bar_width=0.03)
        end
    end
    display(p)
    return p
end