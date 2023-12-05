using Plots, StatsPlots

function plotNNIErrors(mergedists, constraintdiffs, add=false)
    y = Vector{Float64}(mergedists)
    x = Vector{Float64}(sum(constraintdiffs, dims=1)[1,:])
    
    neutral = y .== x
    worse = y .> x
    better = y .< x

    # Jitter the points
    y = y .+ rand(length(y)) / 2
    x = x .+ rand(length(x)) / 2

    if !add
        scatter(x[neutral], y[neutral], xlabel="Sum of errors induced by NNI moves", ylabel="Merged net error", labels="No difference", color="black")
    else
        scatter!(x[neutral], y[neutral], xlabel="Sum of errors induced by NNI moves", ylabel="Merged net error", labels="No difference", color="black")
    end
    scatter!(x[worse], y[worse], labels="Worse than induced", color="red")
    if any(better)
        scatter!(x[better], y[better], labels="Better than induced", color="green")
    end
    plot!([0, maximum(y)], [0, maximum(y)], primary=false)
end

function histNNIErrors(mergedists, constraintdiffs, add=false; kwargs...)
    y = Vector{Float64}(mergedists)
    x = sum(constraintdiffs, dims=1)[1,:]
    y .+= rand(length(y)) / 2

    if !add
        violin(x, y,
            xlabel="Sum of errors induced by NNI moves",
            ylabel="Merged net error"; kwargs...
        )
    else
        violin!(x, y,
            xlabel="Sum of errors induced by NNI moves",
            ylabel="Merged net error"; kwargs...
        )
    end
end


function prettyNNIErrors(dists, constraintdiffs)
    plotNNIErrors(dists, constraintdiffs)
    histNNIErrors(dists, constraintdiffs,
        true, alpha=0.25, color="blue",
        labels="Distributions"
    )
end


function prettyNNIEdgeHeights(dists, constraintdists, edgeheights; metric::Function=median)
    dists = Vector{Float64}(dists)
    constraintdists = Vector{Float64}(sum(constraintdists, dims=1)[1,:])
    edgeheights = Vector{Float64}(metric(edgeheights, dims=1)[1,:])

    x = dists .- constraintdists
    y = edgeheights

    neutral = dists .== constraintdists
    worse = dists .> constraintdists
    better = dists .< constraintdists

    # Jitter the points
    y = y .+ rand(length(y)) / 2
    x = x .+ rand(length(x)) / 5

    scatter(x[neutral], y[neutral], xlabel="Merged network error", ylabel="Metric on heights of NNI move edges", labels="No difference", color="black")
    scatter!(x[worse], y[worse], labels="Worse than induced", color="red")
    scatter!(x[better], y[better], labels="Better than induced", color="green")

    boxplot!((dists .- constraintdists), y, labels=nothing, color="red", alpha=0.25)
end