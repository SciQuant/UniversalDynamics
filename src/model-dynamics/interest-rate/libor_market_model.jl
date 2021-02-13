"""
    LiborMarketModelDynamics{IIP,D,M,DN,T} <: TermStructureModelDynamicsTermStructureModelDynamics{IIP,D,M,DN,T}

Represents the simply-compounded forward interest rates or simple forward rates ``L(t, T,
T+τ)`` under the Libor Market Model.

## Type parameters:
See [`AbstractDynamics`](@ref) for detailed information.

## Fields:
- `attributes`: see [`DynamicsAttributes`](@ref) for detailed information,
- `params`: model parameters, see [`LiborMarketModelParameters`](@ref) for detailed
  information,

## Declaration:

```julia
LiborMarketModelDynamics(
    L0::S, τ, σ, ρ;
    t0=zero(eltype(S)),
    noise::AbstractNoise=DiagonalNoise{length(x0)}(),
    imethod=LiborMarketModelDoNotInterpolate()
    measure::AbstractMeasure=Spot()
) -> LiborMarketModelDynamics
```
"""
struct LiborMarketModelDynamics{IIP,D,M,DN,T,A,P} <: TermStructureModelDynamics{IIP,D,M,DN,T}
    attributes::A
    params::P
end

# see UniversalPricing for τ treatment and implement it here.
function LiborMarketModelDynamics(
    L0::S, τ, σ, ρ;
    noise::AbstractNoise=DiagonalNoise{length(L0)}(),
    measure::AbstractMeasure=Spot(),
    imethod=LiborMarketModelDoNotInterpolate(),
    t0=zero(eltype(S))
) where {S}

    if !(S <: AbstractVector)
        throw(ArgumentError("state *must* be <: AbstractVector."))
    end

    T = eltype(S)
    t0 = convert(T, t0)

    D = length(L0)
    M = dimension(noise)

    # TODO: SPOT rule, check SystemDynamics constructor
    # si es diagonal noise, M tiene que coincidir con D
    if isa(noise, DiagonalNoise) && !isequal(D, M)
        throw(DimensionMismatch("expected `DiagonalNoise` dimension $D, got $M."))
    end

    IIP = isinplace(L0)

    # TODO: SPOT rule, check SystemDynamics constructor
    DN = isa(noise, DiagonalNoise) || (isa(noise, ScalarNoise) && isequal(D, 1))

    # TODO: SPOT rule, check SystemDynamics constructor
    if isequal(ρ, I)
        # should we set or keep ρ to I insted of the following?
        ρ = IIP ? one(T)*I(M) : Diagonal(SVector{M,T}(ones(M)))
    else
        ρsize = (M, M)
        size(ρ) == ρsize || throw(DimensionMismatch("`ρ` *must* be a $(string(ρsize)) matrix."))
        ρ = IIP ? Array{T,2}(ρ) : SMatrix{ρsize...,T}(ρ)
    end

    diffeq_noise = diffeqnoise(t0, ρ, IIP, D, M, DN)
    diffeq_noise_rate_prototype = gprototype(
        LiborMarketModelDynamics{IIP,D,M,DN,T,Nothing,Nothing}(nothing, nothing)
    ) # trick?

    attrs = DynamicsAttributes(t0, L0, ρ, diffeq_noise, diffeq_noise_rate_prototype)

    params = LiborMarketModelParameters{IIP,D,M,DN,T}(t0, L0, τ, σ, ρ, measure, imethod)

    return LiborMarketModelDynamics{IIP,D,M,DN,T,typeof(attrs),typeof(params)}(attrs, params)
end

# IDEA: aca me parece que es conveniente usar traits, ya que tenemos DynamicalSystem,
# ShortRateModelDynamics y SystemDynamics metidos y son todos hijos de AbstractDynamics
for method in (:initialtime, :state, :cor, :noise, :noise_rate_prototype)
    @eval begin
        $method(lmmd::LiborMarketModelDynamics) = $method(lmmd.attributes)
    end
end

# IDEA: usar Traits?
parameters(lmmd::LiborMarketModelDynamics) = lmmd.params

include("libor-market-model/params.jl")
include("libor-market-model/sdes.jl")
