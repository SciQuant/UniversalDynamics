# UniversalDynamics

| **Documentation** |
|:------------ |
| [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://SciQuant.github.io/UniversalDynamics.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://SciQuant.github.io/UniversalDynamics.jl/dev/) |
|**Build Status** |
| [![Build Status](https://github.com/SciQuant/UniversalDynamics.jl/workflows/CI/badge.svg)](https://github.com/SciQuant/UniversalDynamics.jl/actions) |

**UniversalDynamics** is a high-performance library designed to achieve fast and advanced quantitative finance calculations. In few words, it provides useful functionalities for solving Stochastic Differential Equations (SDEs) that are common in quantitative finance. Then, the simulations can be used for pricing financial derivatives with [UniversalPricing.jl](https://github.com/SciQuant/UniversalPricing.jl) solvers.

## Installation

The package can be installed using the Julia package manager. From the Julia REPL, type `]` to enter the `Pkg` REPL mode and run:

```julia
pkg> add https://github.com/SciQuant/UniversalDynamics.jl.git
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add(PackageSpec(url = "https://github.com/SciQuant/UniversalDynamics.jl.git"))
```

## Documentation

- [**STABLE**](https://SciQuant.github.io/UniversalDynamics.jl/stable) &mdash; **Documentation for the most recently tagged version of UniversalDynamics.**
- [**DEVEL**](https://SciQuant.github.io/UniversalDynamics.jl/dev) &mdash; *Documentation for the in-development version of UniversalDynamics.*

## Maintainers and Contributors

**UniversalDynamics** is being developed by [SciQuant](https://github.com/SciQuant), an organization dedicated to creating high-quality scientific software for the financial industry.