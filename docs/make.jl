using Documenter, InPhyNet


makedocs(
    modules=[InPhyNet],
    sitename="InPhyNet.jl",
    authors="Nathan Kolbow",
    clean=true,
    draft=true,
    doctest=false,
    pages=[
        "Introduction" => "index.md",
        "Installation" => "guide/installation.md",
        "Guide" => [
            "How to use InPhyNet" => "guide/how-to-use-inphynet.md",
            "Choosing your subsets" => "guide/subsets.md",
            "Inferring constraint networks" => "guide/constraints.md",
            "Constructing the full network" => "guide/full_net.md"
        ],
        "Simulation Results" => [
            "Runtime evaluations" => "intro/runtimes.md",
            "Accuracy evaluations" => "intro/accuracies.md"
        ],
        "Documentation" => "documentation.md"
    ]
)

deploydocs(
    repo="github.com/NathanKolbow/InPhyNet.jl",
    devurl="stable"
)