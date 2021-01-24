# UniversalDynamics

[![Build Status](https://github.com/SciQuant/UniversalDynamics.jl/workflows/CI/badge.svg)](https://github.com/SciQuant/UniversalDynamics.jl/actions)

UniversalDynamics provides a simple way for defining a Dynamical System formed up by different Dynamics in the financial domain.

```julia
using UniversalDynamics

x0, y0, z0 = rand(1), rand(2), rand(3)
x = SystemDynamics(x0; ρ=ρx, noise=ScalarNoise())
y = SystemDynamics(y0; ρ=ρy) # defaults to DiagonalNoise
z = SystemDynamics(z0; ρ=ρz, noise=NonDiagonalNoise(Mz))

dynamics = (x, y, z)

ds = DynamicalSystem(dynamics)
```