"""
    abstract type VolatilityModelDynamics{IIP,D,M,DN,T} <: ModelDynamics{IIP,D,M,DN,T} end

Supertype for all Volatility Models dynamics, such as:

- [`LocalVolatilityModelDynamics`](@ref),
- [`StochasticVolatilityModelDynamics`](@ref),
"""
abstract type VolatilityModelDynamics{IIP,D,M,DN,T} <: ModelDynamics{IIP,D,M,DN,T} end


struct LocalVolatilityModelDynamics{IIP,D,M,DN,T} <:VolatilityModelDynamics{IIP,D,M,DN,T} end


struct StochasticVolatilityModelDynamics{IIP,D,M,DN,T} <:VolatilityModelDynamics{IIP,D,M,DN,T} end


# should be a subtype of `StochasticVolatilityModel`s?
struct HestonModelDynamics end