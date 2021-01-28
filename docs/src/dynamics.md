## Introduction

In **UniversalDynamics** a *Dynamics* represents continuous time, ``D``-dimensional Ito Systems of Stochastic Differential Equations (SDEs):

```math
d\vec{u}(t) = f(t, \vec{u}(t)) \cdot dt + g(t, \vec{u}(t)) \cdot d\vec{W}(t), \quad \vec{u}(t_0) = \vec{u}_0,\\
```

with drift coefficient ``f \colon \left[t_0, T \right] \times \mathbb{R}^D \rightarrow \mathbb{R}^D``, diffusion coefficient ``g \colon \left[ t_0, T \right] \times \mathbb{R}^D \rightarrow \mathbb{R}^{D \times M}``, ``M``-dimensional driving Wiener correlated or uncorrelated process ``d\vec{W}(t)`` and initial condition ``\vec{u}_0``.

The main abstract type for these kind of objects is given by:

```@docs
UniversalDynamics.AbstractDynamics
```

It is worth mentioning that there exists simpler cases of *Dynamics*, which are fairly common, namely:

```@docs
ScalarNoise
DiagonalNoise
NonDiagonalNoise
```

## Representation

*Dynamics* are represented by two main types, [`SystemDynamics`](@ref) and [`ModelDynamics`](@ref), described below.

### SystemDynamics

```@docs
SystemDynamics
```

### ModelDynamics

```@docs
UniversalDynamics.ModelDynamics
```

[`ModelDynamics`](@ref) subtypes include many common financial models dynamics. Even though they could always be declared as regular [`SystemDynamics`](@ref), it is somewhat useful to have them coded in the library. This enables traceable, reproducible and fast code. Also, for some models, there are many other features implemented in the library, such as Interest Rate Modeling features.

The following is a list with implemented model dynamics:

```@docs
UniversalDynamics.EquityModelDynamics
UniversalDynamics.InterestRateModelDynamics
UniversalDynamics.VolatilityModelDynamics
```

See each corresponding model dynamics section for detailed information.

## Examples

```@example
using UniversalDynamics # hide
using StaticArrays # hide
Dₓ = 3
Mₓ = 5
x0 = @SVector ones(Dₓ)
x = SystemDynamics(x0; noise=NonDiagonalNoise(Mₓ))
```