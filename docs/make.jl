using Pkg
Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter, InPhyNet

makedocs(
    modules=[InPhyNet],
    sitename="InPhyNet.jl",
    authors="Nathan Kolbow",
    clean=true,
    draft=true,
    doctest=false,
    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Walkthrough" => [
            "walkthrough/introduction.md",
            "Estimate a pairwise distance matrix" => "walkthrough/d-matrix.md",
            "Choosing your subsets" => "walkthrough/subsets.md",
            "Inferring constraint networks" => "walkthrough/constraints.md",
            "Constructing the full network" => "walkthrough/full_net.md"
        ],
        # "Simulation Results" => [
        #     "Runtime evaluations" => "intro/runtimes.md",
        #     "Accuracy evaluations" => "intro/accuracies.md"
        # ],
        "Documentation" => "documentation.md"
    ]
)

deploydocs(
    repo="github.com/NathanKolbow/InPhyNet.jl",
    devurl="stable"
)
