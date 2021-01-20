using UniversalDynamics
using Documenter

makedocs(;
    modules=[UniversalDynamics],
    authors="SciQuant",
    repo="https://github.com/rvignolo/UniversalDynamics.jl/blob/{commit}{path}#L{line}",
    sitename="UniversalDynamics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
