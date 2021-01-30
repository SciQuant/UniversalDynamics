## [Introduction](@id ShortRateModelIntroduction)

Under a short rate model, the stochastic state variable is taken to be the instantaneous spot rate. The short rate ``r(t)`` is the (continuously compounded, annualized) interest rate at which an entity can borrow money for an infinitesimally small period of time ``dt`` from time ``t``. Specifying the current short rate does not specify the entire yield curve. However, no-arbitrage arguments show that, under some fairly relaxed technical conditions, if we model the evolution of ``r(t)`` as a stochastic process under the risk-neutral measure ``Q``, then the price at time ``t`` of a zero-coupon or discount bond maturing at time ``T`` with a payoff of 1 is given by:

```math
P(t, T) = E_t^Q \left[ \exp \left( - \int_t^T r(s) ds \right) \right].
```

In principle, the knowledge of the risk neutral dynamics for ``r(t)`` is sufficient to compute the time ``t`` discount bond prices for any maturity ``T > t``. In practice, this expectation may not be computable in closed form, so to make a short rate model operational we must look for subclasses where we can apply fast numerical methods.

**UniversalDynamics** implements the Affine and Quadratic classes, described below.


## Affine Models

```@docs
UniversalDynamics.AffineModelDynamics
```

### One-Factor Affine Model

The spot rate dynamics is given by the following ``1``-dimensional Stochastic Differential Equation in a time span ``\mathbb{I} = \left[ t_0, T \right]``:

```math
dr(t) = \kappa(t) \cdot \left( \theta(t) - r(t) \right) \cdot dt + \Sigma(t) \cdot \sqrt{\alpha(t) + \beta(t) \cdot r(t)} \cdot dW(t),
```

with ``\kappa \colon \mathbb{I} \rightarrow \mathbb{R}``, ``\theta \colon \mathbb{I} \rightarrow \mathbb{R}``, ``\Sigma \colon \mathbb{I} \rightarrow \mathbb{R}``, ``\alpha \colon \mathbb{I} \rightarrow \mathbb{R}``, ``\beta \colon \mathbb{I} \rightarrow \mathbb{R}`` and ``W`` a ``1``-dimensional driving Wiener process in the risk-neutral probability measure ``Q``.

For practical applications of affine models it might be useful to rewrite the model in terms of a factor given by:

```math
r(t) = \xi_0(t) + \xi_1(t) \cdot x(t).
```

with ``\xi_0 \colon \mathbb{I} \rightarrow \mathbb{R}`` and ``\xi_1 \colon \mathbb{I} \rightarrow \mathbb{R}``.

The previous setup yields to zero-coupon bonds of the form:

```math
P(t, T) = \frac{\exp \left( \int_{t_0}^T \xi_0(s) ds \right)}{\exp \left( \int_{t_0}^t \xi_0(s) ds \right)} \cdot \exp \left( A(t, T) - B(t, T) \cdot x(t) \right),
```

where ``A \colon \mathbb{I} × \mathbb{I} → \mathbb{R}`` and ``B \colon \mathbb{I} × \mathbb{I} → \mathbb{R}`` are deterministic functions obtained through the following Riccati System of Ordinary Differential Equations:

```math
\begin{aligned}
  \frac{\partial A}{\partial t} &= \left( \kappa(t) \cdot \theta(t) - \kappa(t) \cdot \xi_0(t) - \xi'_0(t) \right) \cdot \frac{B}{\xi_1(t)} - \frac{1}{2} \cdot \left( \alpha(t) + \beta(t) \cdot \xi_0(t) \right) \cdot \left( \frac{\Sigma(t) \cdot B }{\xi_1(t)} \right)^2, \\
  \frac{\partial B}{\partial t} &= \left( \kappa(t) + \frac{\xi'_1(t)}{\xi_1(t)} \right) \cdot B + \frac{\beta(t)}{2 \cdot \xi_1(t)} \cdot \left( \Sigma(t) \cdot B \right)^2 - \xi_1(t),
\end{aligned}
```

subject to terminal conditions ``A(T, T) = B(T, T) = 0``. Also, notice that when there is no change of variable, i.e. ``\xi_0(t) = 0`` and ``\xi_1(t) = 1``, the discount bond reconstitution formula takes its conventional form:

```math
P(t, T) = \exp \left( A(t, T) - B(t, T) \cdot r(t) \right).
```

On the other hand, if ``\xi_0(t)`` equats to the Instantaneous Forward Rate ``f(t_0, t)`` and ``\xi_1 = 1``, the zero-coupon bond formula becomes:

```math
P(t, T) = \frac{P(t_0, T)}{P(t_0, t)} \cdot \exp \left( A(t, T) - B(t, T) \cdot x(t) \right).
```

Given that ``P(t_0, T)`` corresponds to the initial yield curve, the One-Factor Affine model accepts an additional parameter `P` wich represents ``P(t_0, T) = P_0(T)``. This way, the zero coupon bond evaluation does not need to integrate ``\xi_0``.

Finally, it is worth mentioning that well known models, such as Vasicek, Cox-Ingersoll-Ross, Hull-White and Gaussian Short Rate Models can all be modelled as One-Facto Affine Short Rate Models.

```@docs
UniversalDynamics.OneFactorAffineModelDynamics
```

### Multi-Factor Affine Model

The factor dynamics are given by the following ``D``-dimensional System of Stochastic Differential Equations in a time span ``\mathbb{I} = \left[ t_0, T \right]``:

```math
d\vec{x}(t) = \kappa(t) \cdot \left( \theta(t) - \vec{x}(t) \right) \cdot dt + \Sigma(t) \cdot \mathrm{diag} \left( \sqrt{\alpha(t) + \beta(t) \cdot \vec{x}(t)} \right) \cdot d\vec{W}(t),
```

with ``\kappa \colon \mathbb{I} \rightarrow \mathbb{R}^{D \times D}``, ``\theta \colon \mathbb{I} \rightarrow \mathbb{R}^{D}``, ``\Sigma \colon \mathbb{I} \rightarrow \mathbb{R}^{D \times D}``, ``\alpha \colon \mathbb{I} \rightarrow \mathbb{R}^{D}``, ``\beta \colon \mathbb{I} \rightarrow \mathbb{R}^{D \times D}`` and ``\vec{W}`` a ``D``-dimensional uncorrelated driving Wiener process in the risk-neutral probability measure ``Q``.

The affine short rate is given by:

```math
r(t) = \xi_0(t) + \xi_1(t)^\top \cdot \vec{x}(t).
```

with ``\xi_0 \colon \mathbb{I} \rightarrow \mathbb{R}`` and ``\xi_1 \colon \mathbb{I} \rightarrow \mathbb{R}^D``.

The previous setup yields to zero-coupon bonds of the form:

```math
P(t, T) = \exp \left( A(t, T) - B(t, T)^\top \cdot x(t) \right),
```

where ``A \colon \mathbb{I} × \mathbb{I} → \mathbb{R}`` and ``B \colon \mathbb{I} × \mathbb{I} → \mathbb{R}ᴰ`` are deterministic functions obtained through the following Riccati System of Ordinary Differential Equations:

```math
\begin{aligned}
  \frac{\partial A}{\partial t} &= \theta(t)^\top \cdot \kappa(t)^\top \cdot B - \frac{1}{2} \cdot \alpha(t)^\top \cdot \mathrm{diag} \left( \Sigma(t)^\top \cdot B \right) \cdot \Sigma(t)^\top \cdot B  + \xi_0(t), \\
  \frac{\partial B}{\partial t} &= \kappa(t)^\top \cdot B + \frac{1}{2} \cdot \beta(t)^\top \cdot \mathrm{diag} \left( \Sigma(t)^\top \cdot B \right) \cdot \Sigma(t)^\top \cdot B  + \xi_1(t).
\end{aligned}
```

subject to terminal conditions ``A(T, T) = 0`` and ``B(T, T) = \vec{0}``.

```@docs
UniversalDynamics.MultiFactorAffineModelDynamics
```

## Quadratic Models

```@docs
UniversalDynamics.QuadraticModelDynamics
```

### One-Factor Quadratic Model

The factor dynamics is given by the following ``1``-dimensional Stochastic Differential Equation in a time span ``\mathbb{I} = \left[ t_0, T \right]``:

```math
dx(t) = \kappa(t) \cdot \left( \theta(t) - x(t) \right) \cdot dt + \sigma(t) \cdot dW(t),
```

with ``\kappa \colon \mathbb{I} \rightarrow \mathbb{R}``, ``\theta \colon \mathbb{I} \rightarrow \mathbb{R}``, ``\sigma \colon \mathbb{I} \rightarrow \mathbb{R}`` and ``W`` a ``1``-dimensional driving Wiener process in the risk-neutral probability measure ``Q``.

The quadratic short rate is given by:

```math
r(t) = ξ_0(t) + ξ_1(t) \cdot x(t) + ξ_2(t) \cdot x(t)^2
```

with ``ξ_0 \colon \mathbb{I} \rightarrow \mathbb{R}``, ``ξ_1 \colon \mathbb{I} \rightarrow \mathbb{R}`` and ``ξ_2 \colon \mathbb{I} \rightarrow \mathbb{R}``.

The previous setup yields to zero-coupon bonds of the form:

```math
P(t, T) = \exp \left( -A(t, T) - B(t, T) \cdot x(t) - C(t, T) \cdot x(t)^2 \right),
```

where ``A \colon \mathbb{I} × \mathbb{I} → \mathbb{R}``, ``B \colon \mathbb{I} × \mathbb{I} → \mathbb{R}`` and ``C \colon \mathbb{I} × \mathbb{I} → \mathbb{R}`` are deterministic functions obtained through the following Riccati System of Ordinary Differential Equations:

```math
\begin{aligned}
  \frac{\partial A}{\partial t} &= - \sigma(t)^2 \cdot C - \kappa(t) \cdot \theta(t) \cdot B + \frac{1}{2} \cdot \left( \sigma(t) \cdot B \right)^2 - \xi_0(t), \\
  \frac{\partial B}{\partial t} &= \kappa(t) \cdot B + 2 \cdot \sigma(t)^2 \cdot B \cdot C - 2 \cdot \kappa(t) \cdot \theta(t) \cdot C - \xi_1(t), \\
  \frac{\partial C}{\partial t} &= 2 \cdot \kappa(t) \cdot C + 2 \cdot \left( \sigma(t) \cdot C \right)^2 - \xi_2(t),
\end{aligned}
```

subject to terminal conditions ``A(T, T) = B(T, T) = C(T, T) = 0``.

```@docs
UniversalDynamics.OneFactorQuadraticModelDynamics
```

### Multi-Factor Quadratic Model

The factor dynamics are given by the following ``D``-dimensional System of Stochastic Differential Equations in a time span ``\mathbb{I} = \left[ t_0, T \right]``:

```math
d\vec{x}(t) = \kappa(t) \cdot \left( \theta(t) - \vec{x}(t) \right) \cdot dt + \sigma(t) \cdot d\vec{W}(t),
```

with ``\kappa \colon \mathbb{I} \rightarrow \mathbb{R}^{D \times D}``, ``\theta \colon \mathbb{I} \rightarrow \mathbb{R}^D``, ``\sigma \colon \mathbb{I} \rightarrow \mathbb{R}^{D \times D}`` and ``W`` a ``D``-dimensional uncorrelated driving Wiener process in the risk-neutral probability measure ``Q``.

The quadratic short rate is given by:

```math
r(t) = ξ_0(t) + ξ_1(t)^\top \cdot \vec{x}(t) + \vec{x}(t)^\top \cdot ξ_2(t) \cdot \vec{x}(t)
```

with ``ξ_0 \colon \mathbb{I} \rightarrow \mathbb{R}``, ``ξ_1 \colon \mathbb{I} \rightarrow \mathbb{R}^D`` and ``ξ_2 \colon \mathbb{I} \rightarrow \mathrm{symmetric} \colon \mathbb{R}^{D \times D}``.

The previous setup yields to zero-coupon bonds of the form:

```math
P(t, T) = \exp \left( -A(t, T) - B(t, T)^\top \cdot x(t) - x(t)^\top \cdot C(t, T) \cdot x(t) \right).
```

where ``A \colon \mathbb{I} × \mathbb{I} → \mathbb{R}``, ``B \colon \mathbb{I} × \mathbb{I} → \mathbb{R}^D`` and ``C \colon \mathbb{I} × \mathbb{I} → \mathbb{R}^{D \times D}`` are deterministic functions obtained through the following Riccati System of Ordinary Differential Equations:

```math
\begin{aligned}
  \frac{\partial A}{\partial t} &= - \mathrm{Tr} \left( \sigma(t) \cdot \sigma(t)^\top \cdot C \right) - B^\top \cdot \kappa(t) \cdot \theta(t) + \frac{1}{2} \cdot B^\top \cdot \sigma(t) \cdot \sigma(t)^\top \cdot B - \xi_0(t), \\
  \frac{\partial B}{\partial t} &= \kappa(t)^\top \cdot B + 2 \cdot C \cdot \sigma(t) \cdot \sigma(t)^\top \cdot B - 2 \cdot C \cdot \kappa(t) \cdot \theta(t) - \xi_1(t) \\
  \frac{\partial C}{\partial t} &= \kappa(t)^\top \cdot C + C \cdot \kappa(t) + 2 \cdot C^\top  \cdot \sigma(t) \cdot \sigma(t)^\top \cdot C - \xi_2(t),
\end{aligned}
```

subject to terminal conditions ``A(T, T) = 0``,  ``B(T, T) = \vec{0}`` and  ``C(T, T) = \mathbf{0}``.

```@docs
UniversalDynamics.MultiFactorQuadraticModelDynamics
```
