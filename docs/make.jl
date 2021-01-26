using UniversalDynamics
using Documenter
using DocumenterTools: Themes

# download themes
for file in ("juliadynamics-lightdefs.scss", "juliadynamics-darkdefs.scss", "juliadynamics-style.scss")
    download("https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/$file", joinpath(@__DIR__, file))
end

# create themes
for w in ("light", "dark")
    header = read(joinpath(@__DIR__, "juliadynamics-style.scss"), String)
    theme = read(joinpath(@__DIR__, "juliadynamics-$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "juliadynamics-$(w).scss"), header*"\n"*theme)
end

# compile themes
Themes.compile(joinpath(@__DIR__, "juliadynamics-light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "juliadynamics-dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

makedocs(;
    modules=[UniversalDynamics],
    authors="SciQuant",
    repo="https://github.com/SciQuant/UniversalDynamics.jl/blob/{commit}{path}#L{line}",
    sitename="UniversalDynamics.jl",
    format=Documenter.HTML(;
        prettyurls=false,
        canonical="https://SciQuant.github.io/UniversalDynamics.jl/dev",
        assets = [
        "assets/logo.ico",
        asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
        ],
        collapselevel = 1,
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
