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

@doc raw"""
    MultiFactorAffineModelDynamics{IIP,D,DN,T,S,P,O} <: AffineModelDynamics{MultiFactor,IIP,D,DN,T}

Defines a [`MultiFactor`](@ref) [`ShortRateModel`](@ref) of [`Affine`](@ref) type.
"""
struct MultiFactorAffineModelDynamics{
    IIP, # in or out-of place
    D,   # dimension
    DN,  # has diagonal noise
    T,   # must be eltype(S)
    S,   # state type (initial condition)
    P,   # parameters
    O    # Riccati ODE problem
} <: AffineModelDynamics{MultiFactor,IIP,D,DN,T}
    attributes::A
    params::P
    prob::O
end

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

@doc raw"""
    OneFactorQuadraticModelDynamics{IIP,DN,T,P,O} <: QuadraticModelDynamics{OneFactor,IIP,1,DN,T}

Defines a [`OneFactor`](@ref) [`ShortRateModel`](@ref) of [`Quadratic`](@ref) type.
"""
struct OneFactorQuadraticModelDynamics{IIP,DN,T,P,O} <: QuadraticModelDynamics{OneFactor,IIP,1,DN,T}
    t0::T
    x0::T
    params::P
    prob::O
end

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

@doc raw"""
    MultiFactorQuadraticModelDynamics{IIP,D,DN,T,S,P,O} <: QuadraticModelDynamics{MultiFactor,IIP,D,DN,T}

Defines a [`MultiFactor`](@ref) [`ShortRateModel`](@ref) of [`Quadratic`](@ref) type.
"""
struct MultiFactorQuadraticModelDynamics{
    IIP, # in or out-of place
    D,   # dimension
    DN,  # has diagonal noise
    T,   # must be eltype(S)
    S,   # state type (initial condition)
    P,   # parameters
    O,   # Riccati ODE problem
} <: QuadraticModelDynamics{MultiFactor,IIP,D,DN,T}
    t0::T
    x0::S
    params::P
    prob::O
end

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