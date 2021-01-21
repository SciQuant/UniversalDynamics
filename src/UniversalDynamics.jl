module UniversalDynamics

using StaticArrays
using LinearAlgebra
using DiffEqNoiseProcess

include("auxiliary.jl")

include("noise.jl")
export ScalarNoise, DiagonalNoise, NonDiagonalNoise

include("dynamics.jl")
export SystemDynamics
export isinplace, dimension, noise_dimension, diagonalnoise
export initialtime, state, cor

include("dynamicalsystem.jl")
export DynamicalSystemAttributes , DynamicalSystem

end
