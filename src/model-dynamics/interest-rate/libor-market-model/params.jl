"""
    LiborMarketModelParameters

Wraps all the Libor Market Model relevant attributes in a single object.
"""
struct LiborMarketModelParameters{IIP,D,M,DN,T,
    Q,   # measure
    IT,  # interpolation mode or method
    Te,  # <: AbstractVector
    U,   # <: AbstractVector
    S,   # <: Function that returns an AbstractVector
    R,   # <: AbstractMatrix
} <: InterestRateModelDynamicsParameters
    Tenors::Te
    τ::U
    σ::S
    ρ::R
end
const LMMP = LiborMarketModelParameters

function LiborMarketModelParameters{IIP,D,M,DN,T}(
    t0::Real, L0::AbstractVector, τ::U, σ::S, ρ::R, measure::Q, imethod::IT
) where {IIP,D,M,DN,T,Q,IT,U,S,R}

    Tenors = tenor_structure(τ)
    Te = typeof(Tenors)

    isconcretetype(Q) || error("measure specification error.")
    isconcretetype(IT) || error("interpolation mode specification error.")

    # IDEA: aca se podrian realizar chequeos sobre los argumentos recibidos, por eso dejo
    # disponibles t0 y L0 (para evaluar funciones por ejemplo)

    # caches necesarios tambien podrian ir por aqui, si es que se llega a necesitar alguno

    return LiborMarketModelParameters{IIP,D,M,DN,T,Q,IT,Te,U,S,R}(Tenors, τ, σ, ρ)
end

tenor_structure(τ) = prepend!(cumsum(τ), zero(eltype(τ)))
tenor_structure(τ::SVector) = vcat(SVector((zero(eltype(τ))), cumsum(τ))) # vcat(similar_type(τ, Size(1))(zero(eltype(τ))), cumsum(τ))
