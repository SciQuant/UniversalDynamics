"""
    LiborMarketModelDoNotInterpolate <: LiborMarketModelInterpolationMethod

Avoid interpolation of basic fixed income securities.
"""
struct LiborMarketModelDoNotInterpolate <: LiborMarketModelInterpolationMethod end

interpolate(::BasicFixedIncomeSecurity, ::LiborMarketModelDoNotInterpolate, args...) =
    error("Interpolation of basic fixed income securities is disabled.")