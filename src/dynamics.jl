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


abstract type ModelDynamics{D,M,IIP,DN,T} <: AbstractDynamics{D,M,IIP,DN,T} end


struct SystemDynamics{IIP,D,M,DN,T,S,R} <: AbstractDynamics{IIP,D,M,DN,T}
    t0::T
    x0::S
    ρ::R
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

    return SystemDynamics{IIP,D,M,DN,T,S,R}(t0, x0, ρ)
end

initialtime(sd::SystemDynamics) = sd.t0
state(sd::SystemDynamics) = sd.x0
cor(sd::SystemDynamics{true,D,M,DN,T}) where {D,M,DN,T} = isnothing(sd.ρ) ? one(T)*I(M) : sd.ρ
cor(sd::SystemDynamics{false,D,M,DN,T}) where {D,M,DN,T} = isnothing(sd.ρ) ? Diagonal(SVector{M,T}(ones(M))) : sd.ρ
# ρ(sd::SystemDynamics{false,D,M,DN,T}) where {D,M,DN,T} = isnothing(sd.ρ) ? convert(SMatrix{M,M,T}, I(M)) : sd.ρ
