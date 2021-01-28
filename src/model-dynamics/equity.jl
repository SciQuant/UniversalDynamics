"""
    abstract type EquityModelDynamics{IIP,D,M,DN,T} <: ModelDynamics{IIP,D,M,DN,T} end

Supertype for all Equity Models dynamics, such as:

- [`BlackScholesMertonModelDynamics`](@ref),
-
"""
abstract type EquityModelDynamics{IIP,D,M,DN,T} <: ModelDynamics{IIP,D,M,DN,T} end


"""
    BlackScholesMertonModelDynamics{IIP,D,T} <: EquityModelDynamics{IIP,D,D,true,T}


"""
struct BlackScholesMertonModelDynamics{IIP,D,T} <: EquityModelDynamics{IIP,D,D,true,T} end