
"""
LiborMarketModelInterpolationMethod

Supertype for basic fixed income securities interpolation methods under the Libor Market
Model.
"""
abstract type LiborMarketModelInterpolationMethod end

include("nointerp.jl")
include("schlogl.jl")
# incluir los demas metodos que hay que agregar
