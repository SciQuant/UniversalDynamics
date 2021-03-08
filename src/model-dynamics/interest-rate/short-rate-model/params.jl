"""
    ShortRateParameters{FM,IIP,D,DN} <: InterestRateModelDynamicsParameters

Supertype for [`ShortRateModelDynamics`](@ref) parameters.
"""
abstract type ShortRateParameters{FM,IIP,D,DN} <: InterestRateModelDynamicsParameters end

"""
    AffineParameters{FM,IIP,D,DN} <: ShortRateParameters{FM,IIP,D,DN}

Wraps the model parameters of an [`AffineModelDynamics`](@ref) in a unique object, for
either [`OneFactor`](@ref) and [`MultiFactor`](@ref) cases.
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

function AffineParameters{OneFactor,IIP}(
    t0::Real, r0, κ::K, θ::T, Σ::S, α::A, β::B, ξ₀::X0, ξ₁::X1
) where {IIP,K,T,S,A,B,X0,X1}

    if IIP
        cache = AffineCache(OneFactorAffineModelDynamics{IIP,1,true}, r0)
    else
        cache = nothing
    end
    C = typeof(cache)

    return AffineParameters{OneFactor,IIP,1,true,K,T,S,A,B,X0,X1,C}(κ, θ, Σ, α, β, ξ₀, ξ₁, cache)
end

function AffineParameters{MultiFactor,IIP,D,DN}(
    t0::Real, x0::AbstractVector, κ::K, θ::T, Σ::S, α::A, β::B, ξ₀::X0, ξ₁::X1
) where {IIP,D,DN,K,T,S,A,B,X0,X1}

    if IIP
        cache = AffineCache(MultiFactorAffineModelDynamics{IIP,D,DN}, x0)
    else
        cache = nothing
    end
    C = typeof(cache)

    return AffineParameters{MultiFactor,IIP,D,DN,K,T,S,A,B,X0,X1,C}(κ, θ, Σ, α, β, ξ₀, ξ₁, cache)
end

# OneFactor models do not use caches
function (p::AffineParameters{MultiFactor,true})(t::Real)
    @unpack cache = p
    @unpack κ, θ, Σ, α, β, ξ₁ = cache

    # in place
    p.κ(κ, t)
    p.θ(θ, t)
    p.Σ(Σ, t)
    p.α(α, t)
    p.β(β, t)
    p.ξ₁(ξ₁, t)

    # returns a <: Real
    ξ₀ = p.ξ₀(t)

    return (κ, θ, Σ, α, β, ξ₀, ξ₁)
end

function (p::AffineParameters{FM,false})(t::Real) where {FM}
    return (p.κ(t), p.θ(t), p.Σ(t), p.α(t), p.β(t), p.ξ₀(t), p.ξ₁(t))
end

@doc raw"""
    QuadraticParameters{FM,IIP,D,DN} <: ShortRateParameters{FM,IIP,D,DN}

Wraps the parameters of a [`Quadratic`](@ref) Term Structure (QTS) model in a unique object,
for either [`OneFactor`](@ref) and [`MultiFactor`](@ref) models.
"""
struct QuadraticParameters{FM,IIP,D,DN,K,T,S,X0,X1,X2,C} <: ShortRateParameters{FM,IIP,D,DN}
    κ::K
    θ::T
    σ::S
    ξ₀::X0
    ξ₁::X1
    ξ₂::X2
    cache::C
end

function QuadraticParameters(
    ::Type{OneFactor}, t0::Real, x0, κ::K, θ::T, σ::S, ξ₀::X0, ξ₁::X1, ξ₂::X2
) where {K,T,S,X0,X1,X2}

    IIP = isinplace(x0)
    DN = true
    cache = IIP ? QuadraticCache(OneFactorQuadraticModelDynamics, x0) : nothing
    C = typeof(cache)

    QuadraticParameters{OneFactor,IIP,1,DN,K,T,S,X0,X1,X2,C}(κ, θ, σ, ξ₀, ξ₁, ξ₂, cache)
end

function QuadraticParameters(
    ::Type{MultiFactor}, t0::Real, x0::AbstractVector, κ::K, θ::T, σ::S, ξ₀::X0, ξ₁::X1, ξ₂::X2
) where {K,T,S,X0,X1,X2}

    D = length(x0)
    # pedimos que el usuario explicitamente de un objeto diagonal y no isdiag(σ(t0))
    DN = isa(σ(t0), Diagonal)
    cache = IIP ? QuadraticCache(MultiFactorQuadraticModelDynamics{IIP,D}, x0) : nothing
    C = typeof(cache)

    QuadraticParameters{MultiFactor,IIP,D,DN,K,T,S,X0,X1,X2,C}(κ, θ, σ, ξ₀, ξ₁, ξ₂, cache)
end

(p::QuadraticParameters)(t::Real) = (p.κ(t), p.θ(t), p.σ(t), p.ξ₀(t), p.ξ₁(t), p.ξ₂(t))

struct AffineCache{V,M,VM}
    κ::M
    θ::V
    Σ::VM
    α::V
    β::M
    # ξ₀
    ξ₁::V

    v1::V
    v2::V
    v3::V

    rout::V # riccati out
end

function AffineCache(
    T::Type{<:AffineModelDynamics{FM,IIP,D,DN}}, x::AbstractVector
) where {FM,IIP,D,DN}

    κ = similar(x, D, D)
    θ = similar(x)
    Σ = DN ? similar(x) : similar(x, D, D)
    α = similar(x)
    β = similar(x, D, D)
    ξ₁ = similar(x)

    M = typeof(κ)
    V = typeof(θ)
    VM = typeof(Σ)

    return AffineCache{V,M,VM}(
        κ, θ, Σ, α, β, ξ₁, ntuple(_ -> similar(x), 3)..., similar(x, riccati_dimension(T))
    )
end

struct QuadraticCache{V1,V2,M1,M2,R}
    v1::V1
    v2::V2
    m1::M1
    m2::M2
    rout::R
end

function QuadraticCache(
    T::Type{<:QuadraticModelDynamics{FM,IIP,D}}, x::AbstractVector
) where {FM,IIP,D}
    return MultiFactorQuadraticCache(
        ntuple(_ -> similar(x), 2)...,
        ntuple(_ -> similar(x, D, D), 2)...,
        similar(x, riccati_dimension(T))
    )
end

# factormodel(::ShortRateParameters{FM}) where {FM} = FM
# isinplace(::ShortRateParameters{FM,IIP}) where {FM,IIP} = IIP
# dimension(::ShortRateParameters{FM,IIP,D}) where {FM,IIP,D} = D
diagonalnoise(::ShortRateParameters{FM,IIP,D,DN}) where {FM,IIP,D,DN} = DN