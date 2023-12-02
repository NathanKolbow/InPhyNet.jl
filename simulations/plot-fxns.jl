using Plots

function plotNNIErrors(mergedists, constraintdiffs)
    y = mergedists
    x = sum(constraintdiffs, dims=1)[1,:]
    
    neutral = y .== x
    worse = y .> x
    better = y .< x

    # Jitter the points
    y = y .+ rand(length(dists)) / 2
    x = x .+ rand(length(dists)) / 2

    scatter(x[neutral], y[neutral], xlabel="Sum of induced errors", ylabel="Merged net error", labels="No difference", color="black")
    scatter!(x[worse], y[worse], labels="Worse than induced", color="red")
    if any(better)
        scatter!(x[better], y[better], labels="Better than induced", color="green")
    end
    plot!([0, maximum(y)], [0, maximum(y)], primary=false)
end

