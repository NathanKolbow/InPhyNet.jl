using PhyloNetworks, InPhyNet, Test

test_files = [joinpath(@__DIR__, file) for file in
    [
        "test-AGID-large.jl", "test-AGID.jl", "test-compatibility.jl",
        "test-findvalidpairs.jl", "test-inphynet.jl"
    ]
]

@testset "InPhyNet" begin
    for test_file in test_files
        @testset "$(test_file)" begin
            include(test_file)
        end
    end
end
