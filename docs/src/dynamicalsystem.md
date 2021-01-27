# Dynamical System

Supose we would like to build a System of SDEs that results from the union of a given collection of [`AbstractDynamics`](@ref), with either arbitrary ([`SystemDynamics`](@ref)) or known ([`ModelDynamics`](@ref)) coefficients. To be more precise, given a set of *Dynamics* for ``\{ x(t), y(t), z(t) \}``:

```math
\begin{aligned}
d\vec{x}(t) &= f_x(t, \vec{x}(t)) \cdot dt + g_x(t, \vec{x}(t)) \cdot d\vec{W}_x(t), \quad \vec{x}(t_0) = \vec{x}_0,\\
d\vec{y}(t) &= f_y(t, \vec{y}(t)) \cdot dt + g_y(t, \vec{y}(t)) \cdot d\vec{W}_y(t), \quad \vec{y}(t_0) = \vec{y}_0,\\
d\vec{z}(t) &= f_z(t, \vec{z}(t)) \cdot dt + g_z(t, \vec{z}(t)) \cdot d\vec{W}_z(t), \quad \vec{z}(t_0) = \vec{z}_0,\\
\end{aligned}
```

 we wish to obtain the resulting general *Dynamics* for ``u(t)``:

 ```math
d\vec{u}(t) = f(t, \vec{u}(t)) \cdot dt + g(t, \vec{u}(t)) \cdot d\vec{W}(t), \quad \vec{u}(t_0) = \vec{u}_0,\\
```

with:

```math
\begin{aligned}
\vec{u}(t)       &= \\
f(t, \vec{u}(t)) &= \\
g(t, \vec{u}(t)) &= \\
d\vec{W}(t)      &= \\
\end{aligned}
```



It provides the means to compute solver parameters.

It also provides a binding to many solvers for the resulting System of SDEs.

## Example

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