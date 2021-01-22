module UniversalDynamics

using UnPack: @unpack
using StaticArrays
using LinearAlgebra
using OrdinaryDiffEq
using DiffEqNoiseProcess

include("auxiliary.jl")

include("noise.jl")
export ScalarNoise, DiagonalNoise, NonDiagonalNoise

include("dynamics.jl")
export SystemDynamics
export isinplace, dimension, noise_dimension, diagonalnoise
export initialtime, state, cor

include("securities.jl")
export SystemSecurity

include("dynamicalsystem.jl")
export DynamicalSystem

end
