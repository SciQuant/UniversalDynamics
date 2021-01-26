using UniversalDynamics
using Documenter

makedocs(;
    modules=[UniversalDynamics],
    authors="SciQuant",
    repo="https://github.com/SciQuant/UniversalDynamics.jl/blob/{commit}{path}#L{line}",
    sitename="UniversalDynamics.jl",
    format=Documenter.HTML(;
        prettyurls=false,
        canonical="https://SciQuant.github.io/UniversalDynamics.jl/dev",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "Dynamical systems" => "dynamicalsystem.md"
    ],
)

deploydocs(;
    repo="github.com/SciQuant/UniversalDynamics.jl",
    devbranch = "main"
)
