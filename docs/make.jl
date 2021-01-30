using UniversalDynamics
using Documenter
using DocumenterTools: Themes

# download themes
# for file in ("sciquant-lightdefs.scss", "sciquant-darkdefs.scss", "sciquant-style.scss")
#     download("https://raw.githubusercontent.com/SciQuant/doctheme/master/$file", joinpath(@__DIR__, file))
# end

# create themes
for w in ("light", "dark")
    header = read(joinpath(@__DIR__, "sciquant-style.scss"), String)
    theme = read(joinpath(@__DIR__, "sciquant-$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "sciquant-$(w).scss"), header*"\n"*theme)
end

# compile themes
Themes.compile(joinpath(@__DIR__, "sciquant-light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "sciquant-dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

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
        asset("https://fonts.googleapis.com/css?family=Lato|Source+Code+Pro&display=swap", class=:css),
        ],
        collapselevel = 1,
    ),
    pages=[
        "Introduction" => "index.md",
        "Abstract Dynamics" => [
            "Dynamics" => "ad/dynamics.md",
            "Dynamical system" => "ad/dynamicalsystem.md"
        ],
        "Equity Models" => [
            "Introduction" => "eq/equity.md",
        ],
        "Interest Rate Models" => [
            "Introduction" => "ir/interest_rate.md",
            "Short Rate Models" => "ir/short_rate_model.md",
            "Libor Market Model" => "ir/libor_market_model.md"
        ],
        "Volatility Models" => [
            "Introduction" => "vol/volatility.md",
        ],
        "Simulation" => "simulation.md"
    ],
)

deploydocs(;
    repo="github.com/SciQuant/UniversalDynamics.jl",
    devbranch = "main"
)
