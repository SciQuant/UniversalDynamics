module UniversalDynamics

using UnPack: @unpack
using QuadGK
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
export get_t0, get_state, get_cor, get_noise, get_noise_rate_prototype, get_parameters
export drift, diffusion, drift!, diffusion!

include("securities.jl")
export SystemSecurity
export remake

include("dynamicalsystem.jl")
export DynamicalSystem

include("simulation.jl")
export solve, montecarlo

include("show.jl")

end
