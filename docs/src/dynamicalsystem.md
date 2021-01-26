# Dynamical System

In **UniversalDynamics** a Dynamical System or, more precisely, a Stochastic Dynamical System represents a continuous time, ``D``-dimensional Ito System of Stochastic Differential
Equations:

```math
\begin{aligned}
    d\vec{u}(t)  &= f(t, \vec{u}(t)) \cdot dt + g(t, \vec{u}(t)) \cdot d\vec{W}(t),\\
    \vec{u}(t_0) &= \vec{u}_0,
\end{aligned}
```

where ``f \colon \left[t_0, T \right] \times \mathbb{R}^D \rightarrow \mathbb{R}^D`` represents the drift, ``g \colon \left[ t_0, T \right] \times \mathbb{R}^D \rightarrow \mathbb{R}^{D \times M}`` the diffusion and ``d\vec{W}(t)`` an ``M``-dimensional driving Wiener process.

In the context of quantitative finance we might want to solve a `DynamicalSystem` formed by a set of sub-dynamics, such as:

1. `SystemDynamics` representing arbitrary dynamics;
2. `ModelDynamics` representing known models dynamics, such as:
   - `BlackScholesMerton`,
   - `ShortRateModelDynamics`,
   - `LiborMarketModelDynamics`,
   - `HeathJarrowMortonFrameworkDynamics`,
   - `HestonModelDynamics`,
   - ...

## Example

Supose we want to price a european option on a stock `S` with stochastic interest rates (sacar del cap 1 del Andersen). In this context we need to simulate, for example, a short rate described by any `ShortRateModelDynamics` and a stock price given by a `SystemDynamics`.

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