# import DiffEqBase: isinplace

isinplace(::Type{<:Real}) = false
isinplace(::Real) = false
isinplace(x::AbstractVector) = !isa(x, SVector)
isinplace(x::SubArray{T,1}) where {T} = !isa(x.parent, SVector)