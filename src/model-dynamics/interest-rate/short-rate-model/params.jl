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
    cache = IIP ? AffineCache(OneFactorAffineModelDynamics, r0) : nothing
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

#! Al final es igual al metodo de abajo, asi que unirlos
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