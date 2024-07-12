using Documenter, InPhyNet
using DocThemeIndigo

indigo = DocThemeIndigo.install(InPhyNet)
makedocs(
    modules=[InPhyNet],
    sitename="InPhyNet.jl",
    authors="Nathan Kolbow",
    format=Documenter.HTML(;
        assets=String[indigo]
    ),
    clean=true,
    draft=true,
    doctest=false,
    pages=[
        "Introduction" => "index.md",
        "Installation" => "guide/installation.md",
        "Simulation Study Results" => [
            "Runtime evaluations" => "intro/runtimes.md",
            "Accuracy evaluations" => "intro/accuracies.md"
        ],
        "Guide" => [
            "Using estimated gene trees and SNaQ" => "guide/estgts.md"
        ],
        "Documentation" => "documentation.md"
    ]
)

deploydocs(
    repo="github.com/NathanKolbow/InPhyNet.jl"
)