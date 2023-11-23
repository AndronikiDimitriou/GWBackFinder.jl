using GWBackFinder
using Documenter

DocMeta.setdocmeta!(GWBackFinder, :DocTestSetup, :(using GWBackFinder); recursive=true)

makedocs(;
    modules=[GWBackFinder],
    authors="Androniki Dimitriou",
    repo="https://github.com/AndronikiDimitriou/GWBackFinder.jl/blob/{commit}{path}#{line}",
    sitename="GWBackFinder.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://AndronikiDimitriou.github.io/GWBackFinder.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AndronikiDimitriou/GWBackFinder.jl",
    devbranch="main",
)
