@doc raw"""
    abstract type ShortRateModelDynamics{FM,IIP,D,DN,T} <: TermStructureModelDynamics{IIP,D,D,DN,T} end

Supertype for all Short Rate Models. Notice that the diffusion coefficient for Short Rate
Models is a square matrix ``∈ \mathbb{R}^{D × D}``, i.e. `M = D`.
"""
abstract type ShortRateModelDynamics{FM,IIP,D,DN,T} <: TermStructureModelDynamics{IIP,D,D,DN,T} end

@doc raw"""
    abstract type AffineModelDynamics{FM,IIP,D,DN,T} <: ShortRateModelDynamics{FM,IIP,D,DN,T} end

Supertype for Affine Term Structure (ATS) models. These models provide a closed-form formula
for zero coupon bond prices ``P(t, T)`` in a time span ``\mathbb{I} = \left[ t₀, T
\right]``, given by:
```math
P(t, T) = \exp \left( A(t, T) - B(t, T)^\top \cdot x(t) \right).
```
where ``A \colon \mathbb{I} × \mathbb{I} → \mathbb{R}`` and ``B \colon \mathbb{I} ×
\mathbb{I} → \mathbb{R}ᴰ`` are deterministic functions obtained through a System of Ordinary
Differential Equations called Riccati System.
"""
abstract type AffineModelDynamics{FM,IIP,D,DN,T} <: ShortRateModelDynamics{FM,IIP,D,DN,T} end

@doc raw"""
    abstract type QuadraticModelDynamics{FM,IIP,D,DN,T} <: ShortRateModelDynamics{FM,IIP,D,DN,T} end

Supertype for Quadratic Term Structure (QTS) models. These models provide a closed-form
formula for zero coupon bond prices ``P(t, T)``, given by:
```math
P(t, T) = \exp \left( -A(t, T) - B(t, T)^\top \cdot x(t) - x(t)^\top \cdot C(t, T) \cdot x(t) \right).
```
where ``A(t, T) ``, ``B(t, T)`` and ``C(t, T)`` are deterministic functions obtained through
a System of Ordinary Differential Equations called Riccati System.
"""
abstract type QuadraticModelDynamics{FM,IIP,D,DN,T} <: ShortRateModelDynamics{FM,IIP,D,DN,T} end

@doc raw"""
    FactorModel

Supertype for factor model types.
"""
abstract type FactorModel end

@doc raw"""
    OneFactor <: FactorModel

Encompasses all the [`ShortRateModelDynamics`](@ref)s where a single stochastic factor
``x(t)`` determines the future evolution of all interest rates.
"""
abstract type OneFactor <: FactorModel end

@doc raw"""
    MultiFactor <: FactorModel

Encompasses all the [`ShortRateModelDynamics`](@ref)s where ``N`` stochastic factors
``x(t)`` determine the future evolution of all interest rates.
"""
abstract type MultiFactor <: FactorModel end

factormodel(::ShortRateModelDynamics{FM}) where {FM} = FM

# IDEA: aca me parece que es conveniente usar traits, ya que tenemos DynamicalSystem,
# ShortRateModelDynamics y SystemDynamics metidos y son todos hijos de AbstractDynamics
for method in (:initialtime, :state, :cor, :noise, :noise_rate_prototype)
    @eval begin
        $method(srmd::ShortRateModelDynamics) = $method(srmd.attributes)
    end
end

# IDEA: usar Traits?
parameters(srmd::ShortRateModelDynamics) = srmd.params

include("short-rate-model/models.jl")
include("short-rate-model/params.jl")
include("short-rate-model/sdes.jl")
