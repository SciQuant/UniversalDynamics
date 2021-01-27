# Dynamics

The `AbstractDynamics` type represents continuous time, ``D``-dimensional Ito Systems of Stochastic Differential Equations:

```math
d\vec{u}(t) = f(t, \vec{u}(t)) \cdot dt + g(t, \vec{u}(t)) \cdot d\vec{W}(t), \quad \vec{u}(t_0) = \vec{u}_0,\\
```

with drift coefficient ``f \colon \left[t_0, T \right] \times \mathbb{R}^D \rightarrow \mathbb{R}^D``, diffusion coefficient ``g \colon \left[ t_0, T \right] \times \mathbb{R}^D \rightarrow \mathbb{R}^{D \times M}``, ``M``-dimensional driving Wiener correlated or uncorrelated process ``d\vec{W}(t)`` and initial condition ``\vec{u}_0``.

The previous equation states the most general case of a *Dynamics*, which has non-diagonal noise. There are other simpler noise cases that are really common, namely:

```@docs
ScalarNoise
DiagonalNoise
NonDiagonalNoise
```

## Dynamics representation

*Dynamics* are represented by two main types:

```@docs
SystemDynamics
```

1. `ModelDynamics` representing dynamics with known coefficients.

Dynamics with known coefficients are implemented as subtypes of `ModelDynamics`. Some examples include:

- `BlackScholesMerton`,
- `ShortRateModelDynamics`,
- `LiborMarketModelDynamics`,
- `HeathJarrowMortonFrameworkDynamics`,
- `HestonModelDynamics`,
- ...

# Dynamical System

A `DynamicalSystem` is formed by a collection of `AbstractDynamics`, with either arbitrary or known coefficients.

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