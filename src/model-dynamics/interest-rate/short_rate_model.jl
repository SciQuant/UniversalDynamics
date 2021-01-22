"""
    ShortRateModelDynamics{FM,IIP,D,DN,T} <: TermStructureModelDynamics{IIP,D,D,DN,T}

Supertype for [`ShortRateModel`](@ref)s with [`FactorModel`](@ref) `FM`, dimension `D`,
in-place or out of place version `IIP` and with or without diagonal noise `DS`.
"""
abstract type ShortRateModelDynamics{FM,IIP,D,DN,T} <: TermStructureModelDynamics{IIP,D,D,DN,T} end

@doc raw"""
    AffineModelDynamics{FM,IIP,D,DN,T} <: ShortRateModelDynamics{FM,IIP,D,DN,T}

Supertype for Affine Term Structure (ATS) models. These models provide a closed-form formula
for the zero coupon bond prices `P(t, T)`, given by:
```math
P(t, T) = \exp \left( A(t, T) - B(t, T)^\top \cdot x(t) \right).
```
where ``A(t, T)`` and ``B(t, T)`` are deterministic functions obtained through a System of
Ordinary Differential Equations called Riccati System.
"""
abstract type AffineModelDynamics{FM,IIP,D,DN,T} <: ShortRateModelDynamics{FM,IIP,D,DN,T} end

@doc raw"""
    QuadraticModelDynamics{FM,IIP,D,DN,T} <: ShortRateModelDynamics{FM,IIP,D,DN,T}

Supertype for Quadratic Term Structure (QTS) models. These models provide a closed-form
formula for the zero coupon bond prices `P(t, T)`, given by:
```math
P(t, T) = \exp \left( -A(t, T) - B(t, T)^\top \cdot x(t) - x(t)^\top \cdot C(t, T) \cdot x(t) \right).
```
where ``A(t, T)``, ``B(t, T)`` and ``C(t, T)`` are deterministic functions obtained through
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

Encompasses all the [`ShortRateModel`](@ref)s where a single stochastic factor ``x(t)``
determines the future evolution of all interest rates.
"""
abstract type OneFactor <: FactorModel end

@doc raw"""
    MultiFactor <: FactorModel

Encompasses all the [`ShortRateModel`](@ref)s where ``N`` stochastic factors ``x(t)``
determine the future evolution of all interest rates.
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

parameters(srm::ShortRateModelDynamics) = srm.params # usar Traits?

include("short-rate-model/models.jl")
include("short-rate-model/params.jl")
include("short-rate-model/models.jl")

# include("models.jl")
@doc raw"""
    OneFactorAffineModelDynamics{IIP,T,A,P,O} <: AffineModelDynamics{OneFactor,IIP,1,true,T}

Defines a [`OneFactor`](@ref) [`ShortRateModel`](@ref) of [`Affine`](@ref) type.
"""
struct OneFactorAffineModelDynamics{IIP,T,A,P,O} <: AffineModelDynamics{OneFactor,IIP,1,true,T}
    attributes::A
    params::P
    prob::O
end

function OneFactorAffineModelDynamics(r0::S, κ, θ, Σ, α, β; ξ₀=zero, ξ₁=one, t0=zero(eltype(S))) where {S}

    if !(S <: Union{Real,AbstractVector})
        throw(ArgumentError("state *must* be <: Real/AbstractVector."))
    end

    D = length(r0)
    if !isone(D)
        throw(ArgumentError("state *must* be 1 dimensional."))
    end

    T = eltype(S)
    t0 = convert(T, t0)

    IIP = isinplace(r0)

    params = AffineParameters(OneFactor, t0, r0, κ, θ, Σ, α, β, ξ₀, ξ₁)
    prob = riccati_problem(OneFactorAffineModelDynamics{IIP}, T, params)

    ρ = IIP ? one(T)*I(D) : Diagonal(SVector{D,T}(ones(D)))

    diffeq_noise = diffeqnoise(t0, ρ, IIP, D, D, true) # one factor short rate models have DiagonalNoise
    diffeq_noise_rate_prototype = gprototype(
        OneFactorAffineModelDynamics{IIP,T,Nothing,Nothing,Nothing}(nothing, nothing, nothing)
    ) # trick?

    attrs = DynamicsAttributes(t0, r0, ρ, diffeq_noise, diffeq_noise_rate_prototype)

    A, P, O = typeof.((attrs, params, prob))

    return OneFactorAffineModelDynamics{IIP,T,A,P,O}(attrs, params, prob)
end

# # has r0 field instead of x0
# state(srm::OneFactorAffineModelDynamics) = srm.r0

# @doc raw"""
#     MultiFactorAffineModelDynamics{IIP,D,DN,T,S,P,O} <: AffineModelDynamics{MultiFactor,IIP,D,DN,T}

# Defines a [`MultiFactor`](@ref) [`ShortRateModel`](@ref) of [`Affine`](@ref) type.
# """
# struct MultiFactorAffineModelDynamics{
#     IIP, # in or out-of place
#     D,   # dimension
#     DN,  # has diagonal noise
#     T,   # must be eltype(S)
#     S,   # state type (initial condition)
#     P,   # parameters
#     O    # Riccati ODE problem
# } <: AffineModelDynamics{MultiFactor,IIP,D,DN,T}
#     t0::T
#     x0::S
#     params::P
#     prob::O
# end

# function MultiFactorAffineModelDynamics(x0::S, κ, θ, Σ, α, β, ξ₀, ξ₁; t0=zero(eltype(S))) where {S}

#     if !(S <: AbstractVector)
#         throw(ArgumentError("initial condition *must* be <: AbstractVector."))
#     end

#     T = eltype(S)
#     t0 = convert(T, t0)
#     D = length(x0)
#     IIP = isinplace(x0)
#     params = AffineParameters(MultiFactor, IIP, t0, x0, κ, θ, Σ, α, β, ξ₀, ξ₁)
#     DN = diagonalnoise(params)
#     prob = riccati_problem(MultiFactorAffineModelDynamics{IIP,D}, T, params)

#     P, O = typeof.((params, prob))

#     return MultiFactorAffineModelDynamics{IIP,D,DN,T,S,P,O}(t0, x0, params, prob)
# end

# @doc raw"""
#     OneFactorQuadraticModelDynamics{IIP,DN,T,P,O} <: QuadraticModelDynamics{OneFactor,IIP,1,DN,T}

# Defines a [`OneFactor`](@ref) [`ShortRateModel`](@ref) of [`Quadratic`](@ref) type.
# """
# struct OneFactorQuadraticModelDynamics{IIP,DN,T,P,O} <: QuadraticModelDynamics{OneFactor,IIP,1,DN,T}
#     t0::T
#     x0::T
#     params::P
#     prob::O
# end

# function OneFactorQuadraticModelDynamics(x0::T, κ, θ, σ, ξ₀, ξ₁, ξ₂; t0=zero(S), IIP=false) where {T}

#     if !(T <: Real)
#         throw(ArgumentError("initial condition *must* be <: Real."))
#     end

#     t0 = convert(T, t0)
#     params = QuadraticParameters(OneFactor, IIP, t0, x0, κ, θ, σ, ξ₀, ξ₁, ξ₂)
#     DN = diagonalnoise(params)
#     prob = riccati_problem(OneFactorQuadraticModelDynamics{IIP}, S, params)

#     P, O = typeof.((params, prob))

#     return OneFactorQuadraticModelDynamics{IIP,DN,T,P,O}(t0, x0, params, prob)
# end

# @doc raw"""
#     MultiFactorQuadraticModelDynamics{IIP,D,DN,T,S,P,O} <: QuadraticModelDynamics{MultiFactor,IIP,D,DN,T}

# Defines a [`MultiFactor`](@ref) [`ShortRateModel`](@ref) of [`Quadratic`](@ref) type.
# """
# struct MultiFactorQuadraticModelDynamics{
#     IIP, # in or out-of place
#     D,   # dimension
#     DN,  # has diagonal noise
#     T,   # must be eltype(S)
#     S,   # state type (initial condition)
#     P,   # parameters
#     O,   # Riccati ODE problem
# } <: QuadraticModelDynamics{MultiFactor,IIP,D,DN,T}
#     t0::T
#     x0::S
#     params::P
#     prob::O
# end

# function MultiFactorQuadraticModelDynamics(x0::S, κ, θ, σ, ξ₀, ξ₁, ξ₂; t0=zero(eltype(S))) where {S}

#     if !(S <: AbstractVector)
#         throw(ArgumentError("initial condition *must* be <: AbstractVector."))
#     end

#     T = eltype(S)
#     t0 = convert(T, t0)
#     D = length(x0)
#     IIP = isinplace(x0)
#     params = QuadraticParameters(MultiFactor, IIP, t0, x0, κ, θ, σ, ξ₀, ξ₁, ξ₂)
#     DN = diagonalnoise(params)
#     prob = riccati_problem(MultiFactorQuadraticModelDynamics{IIP,D}, T, params)

#     P, O = typeof.((params, prob))

#     return MultiFactorQuadraticModelDynamics{IIP,D,DN,T,S,P,O}(t0, x0, params, prob)
# end


# include("params.jl")
"""
    ShortRateParameters{FM,D,IIP,DN}

Supertype for [`ShortRateModel`](@ref) parameters with [`FactorModel`](@ref) `FM`, dimension
`D`, in-place or out of place version `IIP` and with or without diagonal noise `DS`.
"""
abstract type ShortRateParameters{FM,IIP,D,DN} <: InterestRateModelDynamicsParameters end

"""
    AffineParameters{FM,D,IIP,DN,K,T,S,A,B,X0,X1,C} <: ShortRateParameters{FM,D,IIP,DN}

Wraps the parameters of an [`Affine`](@ref) Term Structure (ATS) model in a unique object,
for either [`OneFactor`](@ref) and [`MultiFactor`](@ref) models.
"""
struct AffineParameters{FM,IIP,D,DN,K,T,S,A,B,X0,X1,C} <: ShortRateParameters{FM,IIP,D,DN}
    κ::K
    θ::T
    Σ::S
    α::A
    β::B
    ξ₀::X0
    ξ₁::X1
    cache::C
end

function AffineParameters(
    ::Type{OneFactor}, t0::Real, r0, κ::K, θ::T, Σ::S, α::A, β::B, ξ₀::X0, ξ₁::X1
) where {K,T,S,A,B,X0,X1}

    IIP = isinplace(r0)
    DN = true
    cache = IIP ? AffineCache(OneFactorAffineModelDynamics, r0) : nothing # type unstable, para evitarlo podriamos pasar IIP as type?
    C = typeof(cache)

    AffineParameters{OneFactor,IIP,1,DN,K,T,S,A,B,X0,X1,C}(κ, θ, Σ, α, β, ξ₀, ξ₁, cache)
end

function AffineParameters(
    ::Type{MultiFactor}, t0::Real, x0::AbstractVector, κ::K, θ::T, Σ::S, α::A,
    β::B, ξ₀::X0, ξ₁::X1
) where {K,T,S,A,B,X0,X1}

    IIP = isinplace(r0)
    D = length(x0)
    # pedimos que el usuario explicitamente de un objeto diagonal y no isdiag(Σ(t0))
    DN = isa(Σ(t0), Diagonal)
    cache = IIP ? AffineCache(MultiFactorAffineModelDynamics{IIP,D}, x0) : nothing
    C = typeof(cache)

    AffineParameters{MultiFactor,IIP,D,DN,K,T,S,A,B,X0,X1,C}(κ, θ, Σ, α, β, ξ₀, ξ₁, cache)
end

(p::AffineParameters)(t::Real) = (p.κ(t), p.θ(t), p.Σ(t), p.α(t), p.β(t), p.ξ₀(t), p.ξ₁(t))

# @doc raw"""
#     QuadraticParameters{FM,D,IIP,DN,K,T,S,X0,X1,X2,C} <: ShortRateParameters{FM,D,IIP,DN}

# Wraps the parameters of a [`Quadratic`](@ref) Term Structure (QTS) model in a unique object,
# for either [`OneFactor`](@ref) and [`MultiFactor`](@ref) models.
# """
# struct QuadraticParameters{FM,IIP,D,DN,K,T,S,X0,X1,X2,C} <: ShortRateParameters{FM,IIP,D,DN}
#     κ::K
#     θ::T
#     σ::S
#     ξ₀::X0
#     ξ₁::X1
#     ξ₂::X2
#     cache::C
# end

# function QuadraticParameters(
#     ::Type{OneFactor}, IIP::Bool, t0::Real, x0::Real, κ::K, θ::T, σ::S, ξ₀::X0, ξ₁::X1, ξ₂::X2
# ) where {K,T,S,X0,X1,X2}

#     DN = true
#     cache = IIP ? QuadraticCache(OneFactorQuadraticModelDynamics, x0) : nothing
#     C = typeof(cache)

#     QuadraticParameters{OneFactor,IIP,1,DN,K,T,S,X0,X1,X2,C}(κ, θ, σ, ξ₀, ξ₁, ξ₂, cache)
# end

# function QuadraticParameters(
#     ::Type{MultiFactor}, IIP::Bool, t0::Real, x0::AbstractVector, κ::K, θ::T, σ::S, ξ₀::X0,
#     ξ₁::X1, ξ₂::X2
# ) where {K,T,S,X0,X1,X2}

#     D = length(x0)
#     # pedimos que el usuario explicitamente de un objeto diagonal y no isdiag(σ(t0))
#     DN = isa(σ(t0), Diagonal)
#     cache = IIP ? QuadraticCache(MultiFactorQuadraticModelDynamics{IIP,D}, x0) : nothing
#     C = typeof(cache)

#     QuadraticParameters{MultiFactor,IIP,D,DN,K,T,S,X0,X1,X2,C}(κ, θ, σ, ξ₀, ξ₁, ξ₂, cache)
# end

# (p::QuadraticParameters)(t::Real) = (p.κ(t), p.θ(t), p.σ(t), p.ξ₀(t), p.ξ₁(t), p.ξ₂(t))

struct AffineCache{V1,V2,V3,R}
    v1::V1
    v2::V2
    v3::V3
    rout::R
end

function AffineCache(T::Type{<:OneFactorAffineModelDynamics}, x::AbstractVector)
    return AffineCache(
        ntuple(_ -> similar(x), 3)...,
        similar(x, riccati_dimension(T))
    )
end

# function AffineCache(T::Type{<:MultiFactorAffineModelDynamics}, x::AbstractVector)
#     return AffineCache(
#         ntuple(_ -> similar(x), 3)...,
#         similar(x, riccati_dimension(T))
#     )
# end

# struct QuadraticCache{V1,V2,M1,M2,R}
#     v1::V1
#     v2::V2
#     m1::M1
#     m2::M2
#     rout::R
# end

# function QuadraticCache(T::Type{<:OneFactorQuadraticModelDynamics}, ::S) where {S<:Real}
#     return QuadraticCache(
#         ntuple(_ -> nothing, 4)...,
#         similar(Vector{S}, riccati_dimension(T))
#     )
# end

# function QuadraticCache(
#     T::Type{<:MultiFactorQuadraticModelDynamics{IIP,D}}, x::AbstractVector
# ) where {IIP,D}
#     return MultiFactorQuadraticCache(
#         ntuple(_ -> similar(x), 2)...,
#         ntuple(_ -> similar(x, D, D), 2)...,
#         similar(x, riccati_dimension(T))
#     )
# end

# factormodel(::ShortRateParameters{FM}) where {FM} = FM
# isinplace(::ShortRateParameters{FM,IIP}) where {FM,IIP} = IIP
# dimension(::ShortRateParameters{FM,IIP,D}) where {FM,IIP,D} = D
diagonalnoise(::ShortRateParameters{FM,IIP,D,DN}) where {FM,IIP,D,DN} = DN

# include("dynamics.jl")



# esto dummy por ahora aca, va en securities
riccati_dimension(::Type{<:OneFactorAffineModelDynamics}) = 2

function riccati_problem(
    T::Type{<:ShortRateModelDynamics{FM,IIP}}, ::Type{S}, p
) where {FM,IIP,S}
    N = riccati_dimension(T)
    uT = IIP ? zeros(S, N) : SVector{N}(zeros(S, N))
    f = IIP ? riccati! : riccati
    prob =  ODEProblem{IIP}(f, uT, (one(S), zero(S)), p)
    return prob
end

function riccati! end
function riccati end