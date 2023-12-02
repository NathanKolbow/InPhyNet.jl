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