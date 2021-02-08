struct HeathJarrowMortonModelDynamics{IIP,D,M,DN,T,A,P} <: TermStructureModelDynamics{IIP,D,M,DN,T}
    attributes::A
    params::P
end

function HeathJarrowMortonModelDynamics(
    f0::S, τ, σ;
    noise::AbstractNoise=DiagonalNoise{length(f0)}(),
    measure::AbstractMeasure=Spot(),
    t0=zero(eltype(S))
) where {S}

    if !(S <: AbstractVector)
        throw(ArgumentError("state *must* be <: AbstractVector."))
    end

    T = eltype(S)
    t0 = convert(T, t0)

    D = length(f0)
    M = dimension(noise)

    # TODO: SPOT rule, check SystemDynamics constructor
    # si es diagonal noise, M tiene que coincidir con D
    if isa(noise, DiagonalNoise) && !isequal(D, M)
        throw(DimensionMismatch("expected `DiagonalNoise` dimension $D, got $M."))
    end

    IIP = isinplace(f0)

    # TODO: SPOT rule, check SystemDynamics constructor
    DN = isa(noise, DiagonalNoise) || (isa(noise, ScalarNoise) && isequal(D, 1))

    ρ = IIP ? one(T)*I(M) : Diagonal(SVector{M,T}(ones(M)))

    diffeq_noise = diffeqnoise(t0, ρ, IIP, D, M, DN)
    diffeq_noise_rate_prototype = gprototype(
        HeathJarrowMortonModelDynamics{IIP,D,M,DN,T,Nothing,Nothing}(nothing, nothing)
    ) # trick?

    attrs = DynamicsAttributes(t0, f0, ρ, diffeq_noise, diffeq_noise_rate_prototype)

    params = HeathJarrowMortonModelParameters{IIP,D,M,DN,T}(t0, f0, τ, σ, ρ, measure)

    return HeathJarrowMortonModelDynamics{IIP,D,M,DN,T,typeof(attrs),typeof(params)}(attrs, params)
end

# IDEA: aca me parece que es conveniente usar traits, ya que tenemos DynamicalSystem,
# ShortRateModelDynamics y SystemDynamics metidos y son todos hijos de AbstractDynamics
for method in (:initialtime, :state, :cor, :noise, :noise_rate_prototype)
    @eval begin
        $method(hjmd::HeathJarrowMortonModelDynamics) = $method(hjmd.attributes)
    end
end

# IDEA: usar Traits?
parameters(hjmd::HeathJarrowMortonModelDynamics) = hjmd.params

include("heath-jarrow-morton-model/params.jl")
include("heath-jarrow-morton-model/sdes.jl")
