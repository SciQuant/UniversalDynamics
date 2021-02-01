module UniversalDynamics

using UnPack: @unpack
using Memoize
using StaticArrays
using LinearAlgebra
using OrdinaryDiffEq
using StochasticDiffEq
using DiffEqNoiseProcess
using RandomNumbers: Xorshifts

export SVector, SMatrix, @SVector, @SMatrix, Size

include("auxiliary.jl")

include("noise.jl")
export ScalarNoise, DiagonalNoise, NonDiagonalNoise

include("measure.jl")
export Spot, RiskNeutral, Terminal, ForwardTn

include("dynamics.jl")
export SystemDynamics
export isinplace, dimension, noise_dimension, diagonalnoise
export initialtime, state, cor # get_t0, get_state, get_cor

include("securities.jl")
export SystemSecurity
export remake # tendria que expandir el ramake de StochasticDiffEq

include("dynamicalsystem.jl")
export DynamicalSystem

include("simulation.jl")
export solve, montecarlo

include("show.jl")

end
