# Tests for the `inphynet` function in `pipeline.jl`
# function inphynet(
#     estgts::AbstractVector{HybridNetwork},
#     subsets::AbstractVector{<:AbstractVector{<:AbstractString}};
#     verbose::Bool=true,
#     snaqargs...
# )



@testset "Pipeline" begin
    gts = readmultinewick("examples/pipeline-Tgts.tre")
    T = readnewick("examples/pipeline-T.tre")
    subsets = centroid_edge_decomposition(T, 5, 7)

    mnet = inphynet(
        gts,
        subsets
    )
end