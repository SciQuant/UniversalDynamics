import Base: eltype

"""
    abstract type AbstractDynamics{InPlace,Dim,NoiseDim,DiagNoise,elType} end

Supertype for all kind of dynamics.

## Type parameters:
- `InPlace`: states wether coefficients are in or out of place,
- `Dim`: dynamics dimension,
- `NoiseDim`: noise dimension,
- `DiagNoise`: indicates if the noise is of [`DiagonalNoise`](@ref) or
  [`NonDiagonalNoise`](@ref) type,
"""
abstract type AbstractDynamics{InPlace,Dim,NoiseDim,DiagNoise,elType} end

isinplace(::AbstractDynamics{InPlace}) where {InPlace} = InPlace
dimension(::AbstractDynamics{InPlace,Dim}) where {InPlace,Dim} = Dim
noise_dimension(::AbstractDynamics{InPlace,Dim,NoiseDim}) where {InPlace,Dim,NoiseDim} = NoiseDim
diagonalnoise(::AbstractDynamics{InPlace,Dim,NoiseDim,DiagNoise}) where {InPlace,Dim,NoiseDim,DiagNoise} = DiagNoise
eltype(::AbstractDynamics{InPlace,Dim,NoiseDim,DiagNoise,elType}) where {InPlace,Dim,NoiseDim,DiagNoise,elType} = elType

function gprototype(d::AbstractDynamics{IIP,D,M,DN,T}) where {IIP,D,M,DN,T}
    return IIP ? ones(T, gprototype_size(d)...) : @SArray(ones(T, gprototype_size(d)...))
end

gprototype_size(::AbstractDynamics{IIP,1,1}) where {IIP} = (1, ) # either DiagonalNoise or ScalarNoise. Should we return ()? Note that size(::Real) = ()
gprototype_size(::AbstractDynamics{IIP,D,1,false}) where {IIP,D} = (D, ) # ScalarNoise
gprototype_size(::AbstractDynamics{IIP,D,M,true}) where {IIP,D,M} = (D, ) # DiagonalNoise
gprototype_size(::AbstractDynamics{IIP,D,M,false}) where {IIP,D,M} = (D, M) # NonDiagonalNoise


"""
    DynamicsAttributes

Holds parameters related to any [`AbstractDynamics`](@ref).

## Fields:
- `t0`: initial time,
- `x0`: initial state,
- `ρ`: correlation matrix,
- `noise`: Wiener process, and
- `noise_rate_prototype`: prototype or representation of the diffusion coefficient.
"""
struct DynamicsAttributes{T,S,R,N,P}
    t0::T
    x0::S
    ρ::R
    noise::N # diffeq noise always
    noise_rate_prototype::P # it holds my gprototype for SystemDynamics while diffeq gprototype for DynamicalSystem
end

get_t0(attrs::DynamicsAttributes) = attrs.t0
get_state(attrs::DynamicsAttributes) = attrs.x0
get_cor(attrs::DynamicsAttributes) = attrs.ρ
get_noise(attrs::DynamicsAttributes) = attrs.noise
get_noise_rate_prototype(attrs::DynamicsAttributes) = attrs.noise_rate_prototype


"""
    SystemDynamics{IIP,D,M,DN,T} <: AbstractDynamics{IIP,D,M,DN,T}

Represents dynamics with arbitrary coefficients.

## Type parameters:
See [`AbstractDynamics`](@ref) for detailed information.

## Fields:
- `attributes`: see [`DynamicsAttributes`](@ref) for detailed information.

## Declaration

```julia
SystemDynamics(
    x0::S;
    t0=zero(eltype(S)), ρ::R=I, noise::AbstractNoise=DiagonalNoise{length(x0)}(),
) -> SystemDynamics
```

returns a `SystemDynamics` with the given fields, such as state or initial condition `x0`,
intial time `t0`, correlation matrix `ρ` and a driving Wiener process `noise`. Remaining
type parameters are obtained through:

- `IIP`: `true` if `isa(x0, Vector)` or `false` if `isa(x0, Union{Real,SVector})`,
- `D`: equals to `length(x0)`,
- `M`: determined by `noise`, default value is `D`,
- `DN`: determined by `noise`, default value is `true`.
"""
struct SystemDynamics{IIP,D,M,DN,T,A} <: AbstractDynamics{IIP,D,M,DN,T}
    attributes::A
end

function SystemDynamics(
    x0::S; noise::AbstractNoise=DiagonalNoise{length(x0)}(), ρ::R=I, t0=zero(eltype(S))
) where {S,R}

    if !(S <: Union{Real,AbstractVector})
        throw(ArgumentError("state *must* be <: Real/AbstractVector."))
    end

    T = eltype(S)
    t0 = convert(T, t0)

    D = length(x0)
    M = dimension(noise)

    # si es diagonal noise, M tiene que coincidir con D
    if isa(noise, DiagonalNoise) && !isequal(D, M)
        throw(DimensionMismatch("expected `DiagonalNoise` dimension $D, got $M."))
    end

    # Si tenemos un sistema 1D con NonDiagonalNoise, σ es un vector de longitud M. Luego el
    # producto σ ⋅ dW es un vector de longitud 1. Por lo tanto, μ debe ser un vector de
    # longitud 1 y la condicion inicial x0 tambien.
    if S <: Real && isa(noise, NonDiagonalNoise)
        throw(ArgumentError("state *must* be <: AbstractVector for NonDiagonalNoise."))
    end

    IIP = isinplace(x0)

    DN = isa(noise, DiagonalNoise) || (isa(noise, ScalarNoise) && isequal(D, 1))

    if isequal(ρ, I)
        # should we set or keep ρ to I insted of the following?
        ρ = IIP ? one(T)*I(M) : Diagonal(SVector{M,T}(ones(M)))
    else
        ρsize = (M, M)
        size(ρ) == ρsize || throw(DimensionMismatch("`ρ` *must* be a $(string(ρsize)) matrix."))
        ρ = IIP ? Array{T,2}(ρ) : SMatrix{ρsize...,T}(ρ)
    end

    diffeq_noise = diffeqnoise(t0, ρ, IIP, D, M, DN)
    diffeq_noise_rate_prototype = gprototype(SystemDynamics{IIP,D,M,DN,T,Nothing}(nothing)) # trick?

    attrs = DynamicsAttributes(t0, x0, ρ, diffeq_noise, diffeq_noise_rate_prototype)

    return SystemDynamics{IIP,D,M,DN,T,typeof(attrs)}(attrs)
end

for method in (:get_t0, :get_state, :get_cor, :get_noise, :get_noise_rate_prototype)
    @eval begin
        $method(sd::SystemDynamics) = $method(sd.attributes)
    end
end

"""
    abstract type ModelDynamics{IIP,D,M,DN,T} <: AbstractDynamics{IIP,D,M,DN,T} end

Supertype for all dynamics with known coefficients.
"""
abstract type ModelDynamics{IIP,D,M,DN,T} <: AbstractDynamics{IIP,D,M,DN,T} end

include("model-dynamics/equity.jl")
include("model-dynamics/interest_rate.jl")
include("model-dynamics/volatility.jl")
