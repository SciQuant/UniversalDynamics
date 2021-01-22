# UniversalDynamics

[![Build Status](https://github.com/SciQuant/UniversalDynamics.jl/workflows/CI/badge.svg)](https://github.com/SciQuant/UniversalDynamics.jl/actions)

UniversalDynamics provides a simple way for defining a Stochastic Dynamical System formed up by different Stochastic Dynamics.

mmm escribo documentation de una mejor?

```julia
using UniversalDynamics

x0 = rand(1)
x = SystemDynamics(x0; ρ=ρx, noise=ScalarNoise())

y0 = rand(2)
y = SystemDynamics(y0; ρ=ρy) # defaults to DiagonalNoise

z0 = rand(3)
z = SystemDynamics(z0; ρ=ρz, noise=NonDiagonalNoise(Mz))

dynamics = (x, y, z)

ds = DynamicalSystem(dynamics)
```