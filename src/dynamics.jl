import Base: eltype

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


struct DynamicsAttributes{T,S,R,N,P}
    t0::T
    x0::S
    ρ::R
    noise::N # diffeq noise always
    noise_rate_prototype::P # it holds my gprototype for SystemDynamics while diffeq gprototype for DynamicalSystem
end

initialtime(attrs::DynamicsAttributes) = attrs.t0
state(attrs::DynamicsAttributes) = attrs.x0
cor(attrs::DynamicsAttributes) = attrs.ρ
noise(attrs::DynamicsAttributes) = attrs.noise
noise_rate_prototype(attrs::DynamicsAttributes) = attrs.noise_rate_prototype


abstract type ModelDynamics{D,M,IIP,DN,T} <: AbstractDynamics{D,M,IIP,DN,T} end


struct SystemDynamics{IIP,D,M,DN,T,A} <: AbstractDynamics{IIP,D,M,DN,T}
    attributes::A
end

function SystemDynamics(
    x0::S; noise::AbstractNoise=DiagonalNoise{length(x0)}(), ρ::R=nothing, t0=zero(eltype(S))
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

    if isnothing(ρ)
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

initialtime(sd::SystemDynamics) = initialtime(sd.attributes)
state(sd::SystemDynamics) = state(sd.attributes)
cor(sd::SystemDynamics) = cor(sd.attributes)
noise(sd::SystemDynamics) = noise(sd.attributes)
noise_rate_prototype(sd::SystemDynamics) = noise_rate_prototype(sd.attributes)