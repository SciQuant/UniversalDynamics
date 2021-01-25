"""
    LiborMarketModelParameters

Wraps all the Libor Market Model relevant attributes in a single object.
"""
struct LiborMarketModelParameters{IIP,D,M,DN,T,
    E,   # <: AbstractVector
    U,   # <: AbstractVector
    S,   # S might be either <: Function or, if time independent, <: AbstractVector
    R,   # <: AbstractMatrix
    ME,  # measure
    IM   # LMM interpolation method
} <: InterestRateModelDynamicsParameters
    Tenors::E
    τ::U
    σ::S
    ρ::R
    measure::ME
    imethod::IM
end
const LMMP = LiborMarketModelParameters

function LiborMarketModelParameters{IIP,D,M,DN,T}(
    t0::Real, L0::AbstractVector, τ::U, σ::S, ρ::R, measure::ME, imethod::IM
) where {IIP,D,M,DN,T,U,S,R,ME,IM}

    Tenors = tenor_structure(τ)
    E = typeof(Tenors)

    # IDEA: aca se podrian realizar chequeos sobre los argumentos recibidos, por eso dejo
    # disponibles t0 y L0 (para evaluar funciones por ejemplo)

    # caches necesarios tambien podrian ir por aqui, si es que se llega a necesitar alguno

    return LiborMarketModelParameters{IIP,D,M,DN,T,E,U,S,R,ME,MI}(
        Tenors, τ, σ, ρ, measure, imethod
    )
end

tenor_structure(τ) = prepend!(cumsum(τ), zero(eltype(τ)))
tenor_structure(τ::SVector) = vcat(SVector((zero(eltype(τ))), cumsum(τ))) # vcat(similar_type(τ, Size(1))(zero(eltype(τ))), cumsum(τ))
