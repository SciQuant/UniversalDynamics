module UniversalDynamics

using UnPack: @unpack
using Memoize
using StaticArrays
using LinearAlgebra
using OrdinaryDiffEq
using StochasticDiffEq
using DiffEqNoiseProcess

export SVector, SMatrix, @SVector, @SMatrix, Size

include("auxiliary.jl")

include("noise.jl")
export ScalarNoise, DiagonalNoise, NonDiagonalNoise

include("measure.jl")
export Spot, RiskNeutral, Terminal, ForwardTn

include("dynamics.jl")
export SystemDynamics
export isinplace, dimension, noise_dimension, diagonalnoise
export initialtime, state, cor

include("securities.jl")
export SystemSecurity
export remake # tendria que expandir el ramake de StochasticDiffEq

include("dynamicalsystem.jl")
export DynamicalSystem

include("show.jl")

end
