using Plots, StatsPlots

function plotNNIErrors(mergedists, constraintdiffs, add=false; pointalpha=0.75)
    y = Vector{Float64}(mergedists)

    x = Vector{Float64}()
    if size(constraintdiffs, 2) != 1
        x = Vector{Float64}(sum(constraintdiffs, dims=1)[1,:])
    else
        x = Vector{Float64}(constraintdiffs)
    end

    meansx = sort(unique(x))
    means = [mean(y[findall(x .== xval)]) for xval in meansx]
    
    neutral = y .== x
    worse = y .> x
    better = y .< x

    # Jitter the points
    y = y .+ rand(length(y)) / 2
    x = x .+ rand(length(x)) / 2

    if !add
        scatter(x[neutral], y[neutral], xlabel="Sum of errors induced by NNI moves", ylabel="Merged net error", primary=pointalpha > 0., labels="No difference", color="black", alpha=pointalpha)
    else
        scatter!(x[neutral], y[neutral], xlabel="Sum of errors induced by NNI moves", ylabel="Merged net error", primary=pointalpha > 0., labels="No difference", color="black", alpha=pointalpha)
    end
    scatter!(x[worse], y[worse], primary=pointalpha > 0., labels="Worse than induced", color="red", alpha=pointalpha)
    if any(better)
        scatter!(x[better], y[better], primary=pointalpha > 0., labels="Better than induced", color="green", alpha=pointalpha)
    end
    plot!([0, maximum(y)], [0, maximum(y)], primary=false, color="black")
    plot!(meansx, means, labels="Average errors", color="red")
end

function histNNIErrors(mergedists, constraintdiffs, add=false; kwargs...)
    y = Vector{Float64}(mergedists)

    x = Vector{Float64}()
    if size(constraintdiffs, 2) != 1
        x = Vector{Float64}(sum(constraintdiffs, dims=1)[1,:])
    else
        x = Vector{Float64}(constraintdiffs)
    end
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


function prettyNNIErrors(dists, constraintdiffs; pointalpha=0.75)
    plotNNIErrors(dists, constraintdiffs; pointalpha=pointalpha)
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


function plotRobustnessPipelineResults(results::Tuple{Float64, DataFrame, DataFrame})
    allplots = []
    baselineDist, robustDdf, robustNNIdf = results
    
    p = prettyNNIErrors(robustNNIdf[!, "estdists"], robustNNIdf[!, "constraintdists"], pointalpha=0.)
    push!(allplots, p)
    display(p)
    
    b_range = 0:2:(maximum(robustDdf[!,"dists"]))
    colors = ["blue", "green", "red"]
    meanVals = sort(unique(robustDdf[!, "gaussMean"]))
    for (i, meanVal) in enumerate(meanVals)
        fn = histogram!
        if i == 1 fn = histogram end

        p = fn(
            robustDdf[robustDdf[!,"gaussMean"] .== meanVal, "dists"],
            bins=b_range,
            labels="N($(meanVal), $(meanVal))", xlab="Merging error", ylab="Frequency",
            color=colors[i], alpha=0.45,
            normalize=:probability, ylims=[0., 1.]
        )

        if i == 1
            push!(allplots, p)
        end

        if i == length(meanVals)
            display(p)
        end
    end
    return allplots
end