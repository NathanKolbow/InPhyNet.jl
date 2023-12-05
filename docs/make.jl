using Documenter, NetMerge
using DocThemeIndigo

indigo = DocThemeIndigo.install(NetMerge)
makedocs(
    modules=[NetMerge],
    sitename="NetMerge.jl",
    authors="Nathan Kolbow",
    format=Documenter.HTML(;
        assets=String[indigo]
    ),
    clean=true,
    draft=true,
    doctest=false,
    pages=[
        "Introduction" => "index.md",
        "Computational Assessments" => [
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
    repo="github.com/NathanKolbow/network-merging.git",
    devbranch="main"
)