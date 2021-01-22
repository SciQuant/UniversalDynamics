
"""
    InterestRateModelDynamics

Supertype for all Interest Rate Models, such as the [`ShortRateModel`](@ref),
[`LiborMarketModel`](@ref), ForwardMarketModel and Heath-Jarrow-Morton framework.
"""
abstract type InterestRateModelDynamics{IIP,D,M,DN,T} <: ModelDynamics{IIP,D,M,DN,T} end

"""
    TermStructureModelDynamics <: InterestRateModelDynamics

Supertype for all Term Structure Interest Rate models.
"""
abstract type TermStructureModelDynamics{IIP,D,M,DN,T} <: InterestRateModelDynamics{IIP,D,M,DN,T} end

"""
    InterestRateModelDynamicsParameters

Supertype for all Interest Rate Model Dynamics parameters.
"""
abstract type InterestRateModelDynamicsParameters end

include("interest-rate/short_rate_model.jl")
export OneFactorAffineModelDynamics, OneFactorQuadraticModelDynamics,
       MultiFactorAffineModelDynamics, MultiFactorQuadraticModelDynamics

# include("interest-rate/libor_market_model.jl")