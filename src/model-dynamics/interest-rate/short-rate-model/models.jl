@doc raw"""
    OneFactorAffineModelDynamics{IIP,T} <: AffineModelDynamics{OneFactor,IIP,1,true,T}

Represents a [`OneFactor`](@ref) [`AffineModelDynamics`](@ref) type.

## Type parameters:
See [`AbstractDynamics`](@ref) for detailed information.

## Fields:
- `attributes`: see [`DynamicsAttributes`](@ref) for detailed information,
- `params`: model parameters, see [`AffineParameters`](@ref) for detailed information,
- `prob`: Riccati ODEs.

## Declaration:

```julia
OneFactorAffineModelDynamics(
    r0::S, κ, θ, Σ, α, β;
    ξ₀=zero, ξ₁=one, t0=zero(eltype(S))
) -> OneFactorAffineModelDynamics
```

returns a `OneFactorAffineModelDynamics` with the given fields, such as state or initial
condition of the spot rate `r0`, parameters `κ`, `θ`, `Σ`, `α`, `β`, `ξ₀` and `ξ₁` as time
dependent functions and intial time `t0`. `IIP` is the only remaining type parameter left:

- `IIP`: `true` if `isa(x0, Vector)` or `false` if `isa(x0, Union{Real,SVector})`,
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
    MultiFactorAffineModelDynamics{IIP,D,DN,T} <: AffineModelDynamics{MultiFactor,IIP,D,DN,T}

Defines a [`MultiFactor`](@ref) [`AffineModelDynamics`](@ref).

## Type parameters:
See [`AbstractDynamics`](@ref) for detailed information.

## Fields:
- `attributes`: see [`DynamicsAttributes`](@ref) for detailed information,
- `params`: model parameters, see [`AffineParameters`](@ref) for detailed information,
- `prob`: Riccati ODEs.

## Declaration:

```julia
MultiFactorAffineModelDynamics(
    x0::S, κ, θ, Σ, α, β, ξ₀, ξ₁;
    t0=zero(eltype(S))
) -> MultiFactorAffineModelDynamics
```

returns a `MultiFactorAffineModelDynamics` with the given fields, such as state or initial
condition of the factors `x0`, parameters `κ`, `θ`, `Σ`, `α`, `β`, `ξ₀` and `ξ₁` as time
dependent functions and intial time `t0`. Remaining type parameters are obtained through:

- `IIP`: `true` if `isa(x0, Vector)` or `false` if `isa(x0, Union{Real,SVector})`,
- `D`: equals to `length(x0)`,
- `DN`: `true` id `isa(Σ(t0), Diagonal)`.
"""
struct MultiFactorAffineModelDynamics{IIP,D,DN,T,A,P,O} <: AffineModelDynamics{MultiFactor,IIP,D,DN,T}
    attributes::A
    params::P
    prob::O
end

function MultiFactorAffineModelDynamics(x0::S, κ, θ, Σ, α, β, ξ₀, ξ₁; t0=zero(eltype(S))) where {S}

    if !(S <: AbstractVector)
        throw(ArgumentError("state *must* be <: AbstractVector."))
    end

    T = eltype(S)
    t0 = convert(T, t0)

    IIP = isinplace(x0)
    D = length(x0)

    params = AffineParameters(MultiFactor, t0, x0, κ, θ, Σ, α, β, ξ₀, ξ₁)
    prob = riccati_problem(MultiFactorAffineModelDynamics{IIP,D}, T, params)

    DN = diagonalnoise(params)

    ρ = IIP ? one(T)*I(D) : Diagonal(SVector{D,T}(ones(D)))

    diffeq_noise = diffeqnoise(t0, ρ, IIP, D, D, DN)
    diffeq_noise_rate_prototype = gprototype(
        MultiFactorAffineModelDynamics{IIP,D,DN,T,Nothing,Nothing,Nothing}(
            nothing, nothing, nothing
        )
    ) # trick?

    attrs = DynamicsAttributes(t0, x0, ρ, diffeq_noise, diffeq_noise_rate_prototype)

    A, P, O = typeof.((attrs, params, prob))

    return MultiFactorAffineModelDynamics{IIP,D,DN,T,A,P,O}(attrs, params, prob)
end

@doc raw"""
    OneFactorQuadraticModelDynamics{IIP,DN,T,P,O} <: QuadraticModelDynamics{OneFactor,IIP,1,DN,T}

Defines a [`OneFactor`](@ref) [`ShortRateModel`](@ref) of [`Quadratic`](@ref) type.
"""
struct OneFactorQuadraticModelDynamics{IIP,T,A,P,O} <: QuadraticModelDynamics{OneFactor,IIP,1,true,T}
    attributes::A
    params::P
    prob::O
end

function OneFactorQuadraticModelDynamics(x0::S, κ, θ, σ, ξ₀, ξ₁, ξ₂; t0=zero(S)) where {S}

    if !(S <: Union{Real,AbstractVector})
        throw(ArgumentError("state *must* be <: Real/AbstractVector."))
    end

    D = length(r0)
    if !isone(D)
        throw(ArgumentError("state *must* be 1 dimensional."))
    end

    T = eltype(S)
    t0 = convert(T, t0)

    IIP = isinplace(x0)

    params = QuadraticParameters(OneFactor, t0, x0, κ, θ, σ, ξ₀, ξ₁, ξ₂)
    prob = riccati_problem(OneFactorQuadraticModelDynamics{IIP}, T, params)

    ρ = IIP ? one(T)*I(D) : Diagonal(SVector{D,T}(ones(D)))

    diffeq_noise = diffeqnoise(t0, ρ, IIP, D, D, true) # one factor short rate models have DiagonalNoise
    diffeq_noise_rate_prototype = gprototype(
        OneFactorQuadraticModelDynamics{IIP,T,Nothing,Nothing,Nothing}(nothing, nothing, nothing)
    ) # trick?

    attrs = DynamicsAttributes(t0, r0, ρ, diffeq_noise, diffeq_noise_rate_prototype)

    A, P, O = typeof.((attrs, params, prob))

    return OneFactorQuadraticModelDynamics{IIP,T,A,P,O}(attrs, params, prob)
end

@doc raw"""
    MultiFactorQuadraticModelDynamics{IIP,D,DN,T,S,P,O} <: QuadraticModelDynamics{MultiFactor,IIP,D,DN,T}

Defines a [`MultiFactor`](@ref) [`ShortRateModel`](@ref) of [`Quadratic`](@ref) type.
"""
struct MultiFactorQuadraticModelDynamics{IIP,D,DN,T,A,P,O,} <: QuadraticModelDynamics{MultiFactor,IIP,D,DN,T}
    attributes::A
    params::P
    prob::O
end

function MultiFactorQuadraticModelDynamics(x0::S, κ, θ, σ, ξ₀, ξ₁, ξ₂; t0=zero(eltype(S))) where {S}

    if !(S <: AbstractVector)
        throw(ArgumentError("state *must* be <: AbstractVector."))
    end

    T = eltype(S)
    t0 = convert(T, t0)

    IIP = isinplace(x0)
    D = length(x0)

    params = QuadraticParameters(MultiFactor, t0, x0, κ, θ, σ, ξ₀, ξ₁, ξ₂)
    prob = riccati_problem(MultiFactorQuadraticModelDynamics{IIP,D}, T, params)

    DN = diagonalnoise(params)

    ρ = IIP ? one(T)*I(D) : Diagonal(SVector{D,T}(ones(D)))

    diffeq_noise = diffeqnoise(t0, ρ, IIP, D, D, DN)
    diffeq_noise_rate_prototype = gprototype(
        MultiFactorQuadraticModelDynamics{IIP,D,DN,T,Nothing,Nothing,Nothing}(
            nothing, nothing, nothing
        )
    ) # trick?

    attrs = DynamicsAttributes(t0, x0, ρ, diffeq_noise, diffeq_noise_rate_prototype)

    A, P, O = typeof.((attrs, params, prob))

    return MultiFactorQuadraticModelDynamics{IIP,D,DN,T,A,P,O}(attrs, params, prob)
end