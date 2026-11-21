# To build locally, from the repo root:
#   julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#   julia --project=docs docs/make.jl
using Documenter, InPhyNet

makedocs(
    modules=[InPhyNet],
    sitename="InPhyNet.jl",
    authors="Nathan Kolbow",
    format=Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link = "main",
    ),
    clean=true,
    doctest=false,
    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "How InPhyNet works" => "method.md",
        "Walkthrough" => [
            "walkthrough/introduction.md",
            "Estimate a pairwise distance matrix" => "walkthrough/d-matrix.md",
            "Choosing your subsets" => "walkthrough/subsets.md",
            "Inferring constraint networks" => "walkthrough/constraints.md",
            "Constructing the full network" => "walkthrough/full_net.md"
        ],
        "Empirical example" => "empirical.md",
        "Advanced" => [
            "Pairwise InPhyNet" => "advanced/pairwise-usage.md"
        ],
        "Documentation" => "documentation.md"
    ]
)

deploydocs(;
    repo="github.com/NathanKolbow/InPhyNet.jl",
    devurl="stable",
    devbranch="main",
)
