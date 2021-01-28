## Introduction

Supose we would like to build a System of SDEs that comes from the union of a given collection of [`AbstractDynamics`](@ref), with either arbitrary (see [`SystemDynamics`](@ref)) or known (see [`ModelDynamics`](@ref)) coefficients. To be more precise, given a set of *Dynamics* for ``\{ x(t), y(t), z(t) \}``:

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
\vec{u}(t) =
    \begin{bmatrix}
        \vec{x}(t) \\
        \vec{y}(t) \\
        \vec{z}(t)
    \end{bmatrix}
\quad
f(t, \vec{u}(t)) =
    \begin{bmatrix}
        f_x(t, \vec{x}(t)) \\
        f_y(t, \vec{y}(t)) \\
        f_z(t, \vec{z}(t))
    \end{bmatrix}
\quad
g(t, \vec{u}(t)) =
    \begin{bmatrix}
        g_x(t, \vec{x}(t)) & 0 & 0 \\
        0 & g_y(t, \vec{y}(t)) & 0 \\
        0 & 0 & g_z(t, \vec{z}(t))
    \end{bmatrix}
\quad
d\vec{W}(t) =
    \begin{bmatrix}
        \vec{W}_x(t) \\
        \vec{W}_y(t) \\
        \vec{W}_z(t)
    \end{bmatrix}
```

A [`DynamicalSystem`](@ref) provides a shorthand for constructing all the previous prototypes which are needed in the Stochastic Differential Equation solvers. However, it does not construct the general functions ``f`` and ``g``. This task is left to the user. It does provide many useful handlers that come in handy when coding ``f`` and ``g``.

```@docs
DynamicalSystem
```

## Simple example

```@example
using OrderedCollections # hide
using UniversalDynamics # hide
# declare dynamics
x = SystemDynamics(rand(1); noise=ScalarNoise())
y = SystemDynamics(rand(2); noise=DiagonalNoise(2), ρ=[1 0.3; 0.3 1])
z = SystemDynamics(rand(3); noise=NonDiagonalNoise(2), ρ=[1 0.2; 0.2 1])

# group dynamics in a container
dynamics = OrderedDict(:x => x, :y => y, :z => z)

# compute dynamical system
ds = DynamicalSystem(dynamics)
```

## Out of place example

```@example
using OrderedCollections # hide
using UniversalDynamics # hide
using StaticArrays # hide
using UnPack # hide
# load some parameters
include("../../test/DaiSingletonParameters_A3_1.jl")

# define short rate model dynamics parameters
x0 = @SVector [υ₀, θ₀, r₀]

ξ₀(t) = zero(t) # ξ₀ = zero

ξ₁(t) = @SVector [0, 0, 1]

ϰ(t) = @SMatrix([
    μ     0 0
    0     ν 0
    κ_rυ -κ κ
])

θ(t) = @SVector [ῡ, θ̄, θ̄ ]

Σ(t) = @SMatrix [
    η           0    0
    η * σ_θυ    1 σ_θr
    η * σ_rυ σ_rθ    1
]

α(t) = @SVector [0, ζ^2, α_r]

β(t) = @SMatrix [
    1   0 0
    β_θ 0 0
    1   0 0
]

# declare short rate model dynamics
x = MultiFactorAffineModelDynamics(x0, ϰ, θ, Σ, α, β, ξ₀, ξ₁)

# declare money market account dynamics
B = SystemDynamics(one(eltype(x)))

# out of place drift coefficient
function f(u, p, t)
    @unpack x_dynamics, x_security, B_security = p

    x = remake(x_security, u)
    B = remake(B_security, u)

    IR = FixedIncomeSecurities(x_dynamics, x, B)

    dx = drift(x(t), parameters(x_dynamics), t)
    dB = IR.r(t) * B(t)

    return vcat(dx, dB)
end

# out of place diffusion coefficient
function g(u, p, t)
    @unpack x_dynamics, x_security, B_security = p

    x = remake(x_security, u)
    B = remake(B_security, u)

    dx = diffusion(x(t), parameters(x_dynamics), t)
    dB = zero(eltype(u))

    return @SMatrix [dx[1,1] dx[1,2] dx[1,3]  0
                     dx[2,1] dx[2,2] dx[2,3]  0
                     dx[3,1] dx[3,2] dx[3,3]  0
                           0       0       0 dB]
end

# group dynamics in a container
dynamics = OrderedDict(:x => x, :B => B)

# compute dynamical system
ds = DynamicalSystem(f, g, dynamics, nothing)
```
