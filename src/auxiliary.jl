# import DiffEqBase: isinplace

isinplace(::Type{<:Real}) = false
isinplace(::Real) = false
isinplace(x::AbstractVector) = !isa(x, SVector)
isinplace(x::SubArray) = isinplace(x.parent)
