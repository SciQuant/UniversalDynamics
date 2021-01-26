# Dynamical System

In **UniversalDynamics** a Dynamical System or, more precisely, a Stochastic Dynamical System represents a continuous time, ``D``-dimensional Ito System of Stochastic Differential
Equations:

```math
d\vec{u}(t) = f(t, \vec{u}(t)) \cdot dt + g(t, \vec{u}(t)) \cdot d\vec{W}(t), \quad \vec{u}(t_0) = \vec{u}_0,\\
```

with drift coefficient ``f \colon \left[t_0, T \right] \times \mathbb{R}^D \rightarrow \mathbb{R}^D``, diffusion coefficient ``g \colon \left[ t_0, T \right] \times \mathbb{R}^D \rightarrow \mathbb{R}^{D \times M}``, ``M``-dimensional driving Wiener correlated or uncorrelated process ``d\vec{W}(t)`` and initial condition ``\vec{u}_0``.

The previous equation represent the most general case of a Dynamical System, which is referenced as the non-diagonal noise case. There are other simpler cases that are really important and are implemented in the library, namely:

```@docs
ScalarNoise
DiagonalNoise
```

In the context of quantitative finance we might want to declare a `DynamicalSystem` formed by a set of dynamics, with either arbitrary or known coefficients. To take this into account, **UniversalDynamics** uses the following type architecture:

```julia
abstract type AbstractDynamics{InPlace,Dim,NoiseDim,DiagNoise,elType} end

abstract type ModelDynamics{D,M,IIP,DN,T} <: AbstractDynamics{D,M,IIP,DN,T} end

struct SystemDynamics{IIP,D,M,DN,T,A} <: AbstractDynamics{IIP,D,M,DN,T}
    attributes::A
end
```

```@docs
SystemDynamics
```

1. `SystemDynamics` representing arbitrary dynamics;
2. `ModelDynamics` representing known models dynamics, such as:
   - `BlackScholesMerton`,
   - `ShortRateModelDynamics`,
   - `LiborMarketModelDynamics`,
   - `HeathJarrowMortonFrameworkDynamics`,
   - `HestonModelDynamics`,
   - ...

## Dynamical System definition

Supose we want to price a european option on a stock `S` with stochastic interest rates (sacar del cap 1 del Andersen). In this context we need to simulate, for example, a short rate described by any `ShortRateModelDynamics` and a stock price given by a `SystemDynamics`.

TODO: equations

```julia
# define dynamics
x = MultiFactorAffineModelDynamics(x0, ϰ, θ, Σ, α, β, ξ₀, ξ₁)
S = SystemDynamics(S0; noise=NonDiagonalNoise(Mₛ))

# container
dynamics = OrderedDict(:x => x, :S => S)

# define dynamical system formed by the given dynamics
dynamical_system = DynamicalSystem(dynamics)
```

This will allow the user to check important attributes...

However, in order to solve a `DynamicalSystem`, the drift `f` and the diffusion `g` functions must be provided.