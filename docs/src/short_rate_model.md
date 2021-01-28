## Introduction

Under a short rate model, the stochastic state variable is taken to be the instantaneous spot rate. The short rate ``r(t)`` is the (continuously compounded, annualized) interest rate at which an entity can borrow money for an infinitesimally small period of time ``dt`` from time ``t``. Specifying the current short rate does not specify the entire yield curve. However, no-arbitrage arguments show that, under some fairly relaxed technical conditions, if we model the evolution of ``r(t)`` as a stochastic process under the risk-neutral measure ``Q``, then the price at time ``t`` of a zero-coupon or discount bond maturing at time ``T`` with a payoff of 1 is given by:
```math
P(t, T) = E_t^Q \left[ \exp \left( - \int_t^T r(s) ds \right) \right].
```

In principle, the knowledge of the risk neutral dynamics for ``r(t)`` is sufficient to compute the time ``t`` discount bond prices for any maturity ``T > t``. In practice, this expectation may not be computable in closed form, so to make a short rate model operational we must look for subclasses where we can apply fast numerical methods.

The main abstract type for short rate model dynamics is given by:

```@docs
UniversalDynamics.ShortRateModelDynamics
```

**UniversalDynamics** implements the Affine and Quadratic clases, described below.


## Affine Models

```@docs
UniversalDynamics.AffineModelDynamics
```

### One-Factor Affine Model

ir al andersen y ver

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

The previous setup yields to zero-coupon bonds:

```math
P(t, T) = \exp \left( A(t, T) - B(t, T)^\top \cdot x(t) \right),
```

where ``A \colon \mathbb{I} × T → \mathbb{R}`` and ``B \colon \mathbb{I} × T → \mathbb{R}ᴰ``
are deterministic functions obtained through the following Riccati System of Ordinary Differential Equations:

```math
\begin{aligned}
  \frac{\partial A}{\partial t} &= \theta(t)^\top \cdot \kappa(t)^\top \cdot B - \frac{1}{2} \cdot \alpha(t)^\top \cdot \mathrm{diag} \left( \Sigma(t)^\top \cdot B \right) \cdot \Sigma(t)^\top \cdot B  + \xi_0(t), \\
  \frac{\partial B}{\partial t} &= \kappa(t)^\top \cdot B + \frac{1}{2} \cdot \beta(t)^\top \cdot \mathrm{diag} \left( \Sigma(t)^\top \cdot B \right) \cdot \Sigma(t)^\top \cdot B  + \xi_1(t).
\end{aligned}
```

subject to terminal conditions ``A(T, T) = 0`` and ``B(T, T) = 0``.


```@docs
UniversalDynamics.MultiFactorAffineModelDynamics
```

## Quadratic Models

```@docs
UniversalDynamics.QuadraticModelDynamics
```