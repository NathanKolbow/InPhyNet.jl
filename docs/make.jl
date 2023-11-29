using Documenter, NetMerge
using DocThemeIndigo

indigo = DocThemeIndigo.install(NetMerge)
makedocs(
    modules=[NetMerge],
    sitename="NetMerge.jl",
    authors="Nathan Kolbow",
    format=Documenter.HTML(;
        prettyurls=true,
        assets=String[indigo]
    ),
    clean=true,
    repo="NetMerge",
    draft=true,
    pages=[
        "Introduction" => "index.md",
        "Guide" => [
            "Using estimated gene trees" => "guide.md"
        ],
        "Documentation" => "documentation.md"
    ],
    checkdocs=:exports
)

deploydocs(
    repo="github.com/NathanKolbow/network-merging.git"
)