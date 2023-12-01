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

# Waiting for the repo to be public to do deployment
#
# deploydocs(
#     devbranch="gh-pages",
#     repo="github.com/NathanKolbow/network-merging.git"
# )