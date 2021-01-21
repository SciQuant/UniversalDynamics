using UniversalDynamics
using Test

using DiffEqNoiseProcess
using LinearAlgebra
using StaticArrays

@testset "Out of Place" begin

    @testset "OneDimensional + ScalarNoise" begin
        IIP = false
        D = 1
        M = 1
        DN = true
        T = Float64

        t0 = zero(T)
        x0 = T(0.1) # SVector(0.1)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=ScalarNoise())

        @test isinplace(x_dynamics) == IIP
        @test dimension(x_dynamics) == D
        @test noise_dimension(x_dynamics) == M
        @test diagonalnoise(x_dynamics) == DN
        @test UniversalDynamics.eltype(x_dynamics) == T
        @test initialtime(x_dynamics) == t0
        @test state(x_dynamics) == x0
        @test cor(x_dynamics) == Diagonal(SVector{M,T}(ones(M))) == I

        attrs = DynamicalSystemAttributes((x_dynamics, ))
        @test initialtime(attrs) == t0
        @test state(attrs) == SVector(x0) # once we define a system, we use vectors. Check!
        @test cor(attrs) == I
        @test isnothing(attrs.noise)
        @test isnothing(attrs.noise_rate_prototype) # define gprototype?

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        @test isinplace(ds) == IIP
        @test dimension(ds) == D
        @test noise_dimension(ds) == M
        @test diagonalnoise(ds) == DN
        @test UniversalDynamics.eltype(ds) == T
    end

    # should yield to same results as the "OneDimensional + ScalarNoise" test
    @testset "OneDimensional + DiagonalNoise" begin
        IIP = false
        D = 1
        M = 1
        DN = true
        T = Float64

        t0 = zero(T)
        x0 = T(0.1) # SVector(0.1)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=DiagonalNoise(M))

        @test isinplace(x_dynamics) == IIP
        @test dimension(x_dynamics) == D
        @test noise_dimension(x_dynamics) == M
        @test diagonalnoise(x_dynamics) == DN
        @test UniversalDynamics.eltype(x_dynamics) == T
        @test initialtime(x_dynamics) == t0
        @test state(x_dynamics) == x0
        @test cor(x_dynamics) == Diagonal(SVector{M,T}(ones(M))) == I

        attrs = DynamicalSystemAttributes((x_dynamics, ))
        @test initialtime(attrs) == t0
        @test state(attrs) == SVector(x0) # once we define a system, we use vectors. Check!
        @test cor(attrs) == I
        @test isnothing(attrs.noise)
        @test isnothing(attrs.noise_rate_prototype) # define gprototype?

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        @test isinplace(ds) == IIP
        @test dimension(ds) == D
        @test noise_dimension(ds) == M
        @test diagonalnoise(ds) == DN
        @test UniversalDynamics.eltype(ds) == T
    end

    @testset "OneDimensional + NonDiagonalNoise" begin
        IIP = false
        D = 1
        M = 3
        DN = false
        T = Float64

        t0 = zero(T)
        x0 = 0.1
        @test_throws ArgumentError SystemDynamics(x0; t0=t0, noise=NonDiagonalNoise(M))
        x0 = SVector{1,T}(0.1) # must be a Vector
        x_dynamics = SystemDynamics(x0; t0=t0, noise=NonDiagonalNoise(M))

        @test isinplace(x_dynamics) == IIP
        @test dimension(x_dynamics) == D
        @test noise_dimension(x_dynamics) == M
        @test diagonalnoise(x_dynamics) == DN
        @test UniversalDynamics.eltype(x_dynamics) == T
        @test initialtime(x_dynamics) == t0
        @test state(x_dynamics) == x0
        @test cor(x_dynamics) == Diagonal(SVector{M,T}(ones(M))) == I

        attrs = DynamicalSystemAttributes((x_dynamics, ))
        @test initialtime(attrs) == t0
        @test state(attrs) == x0
        @test cor(attrs) == I
        @test attrs.noise == WienerProcess(attrs.t0, @SVector(zeros(M)), @SVector(zeros(M)))
        @test attrs.noise_rate_prototype == @SMatrix(ones(D, M))

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        @test isinplace(ds) == IIP
        @test dimension(ds) == D
        @test noise_dimension(ds) == M
        @test diagonalnoise(ds) == DN
        @test UniversalDynamics.eltype(ds) == T
    end


    @testset "MultiDimensional + ScalarNoise" begin
        IIP = false
        D = 3
        M = 1
        DN = false
        T = Float64

        t0 = zero(T)
        x0 = @SVector rand(D)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=ScalarNoise())

        @test isinplace(x_dynamics) == IIP
        @test dimension(x_dynamics) == D
        @test noise_dimension(x_dynamics) == M
        @test diagonalnoise(x_dynamics) == DN
        @test UniversalDynamics.eltype(x_dynamics) == T
        @test initialtime(x_dynamics) == t0
        @test state(x_dynamics) == x0
        @test cor(x_dynamics) == Diagonal(SVector{M,T}(ones(M))) == I

        attrs = DynamicalSystemAttributes((x_dynamics, ))
        @test initialtime(attrs) == t0
        @test state(attrs) == x0
        @test cor(attrs) == I
        @test attrs.noise == WienerProcess(attrs.t0, @SVector(zeros(1)), @SVector(zeros(1)))
        @test isnothing(attrs.noise_rate_prototype) # diffeq will use u0 as prototype!

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        @test isinplace(ds) == IIP
        @test dimension(ds) == D
        @test noise_dimension(ds) == M
        @test diagonalnoise(ds) == DN
        @test UniversalDynamics.eltype(ds) == T
    end
end

# llamar a los test Mixed Dimensions o algo asi cuando empiece a combinar


# Case of study:
#  - OneDimensional -> D == 1
#  - ScalarNoise or DiagonalNoise -> M == D == 1
#  - Out Of Place -> IIP = false
#  - Without Interest Rates
# SDEs:
#   dx(t) = r_x ⋅ x(t) ⋅ dt + σ_x ⋅ x(t) ⋅ dW_x(t)

# x is a one dimensional security process with DiagonalNoise
x0 = 0.1 # SVector(0.1)
x_dynamics = SystemDynamics(x0; t0=zero(x0), noise=DiagonalNoise(1))

attrs = DynamicalSystemAttributes((x_dynamics, ))

@test attrs.x0 == SVector(x0)
@test attrs.ρ == I
@test isnothing(attrs.noise)
@test isnothing(attrs.noise_rate_prototype)

# Case of study:
#   Same as above but using ScalarNoise() instead of DiagonalNoise(1)

# x is a one dimensional security process with ScalarNoise (i.e. equivalent to DiagonalNoise)
x0 = 0.1 # SVector(0.1)
x_dynamics = SystemDynamics(x0; t0=zero(x0), m=ScalarNoise())

attrs = DynamicalSystemAttributes((x_dynamics, ))

@test attrs.x0 == SVector(x0)
@test attrs.ρ == I
@test isnothing(attrs.noise)
@test isnothing(attrs.noise_rate_prototype)

# Case of study:
#  - OneDimensional -> D == 1
#  - NonDiagonalNoise with M > D
#  - Out Of Place -> IIP = false
#  - Without Interest Rates
# SDEs:
#   dx(t) = μ_x(x(t), t) ⋅ dt + σ_x(x(t), t) ⋅ dW_x(t)

# x is a one dimensional security process with NonDiagonalNoise (i.e. equivalent to DiagonalNoise)
D = 1
M = 3
x0 = SVector(0.1) # no puede ser <: Real ya que si o si el producto σ_x * dW_x retorna un Vector
x_dynamics = SystemDynamics(x0; t0=zero(eltype(x0)), m=NonDiagonalNoise(M))

attrs = DynamicalSystemAttributes((x_dynamics, ))

@test isinplace(attrs) == false
@test dimension(attrs) == D
@test noise_dimension(attrs) == M
@test diagonalnoise(attrs) == false
@test attrs.x0 == SVector(x0)
@test attrs.ρ == I
@test attrs.noise == WienerProcess(attrs.t0, @SVector(zeros(3)), @SVector(zeros(3)))
@test attrs.noise_rate_prototype == @SMatrix(ones(1, 3))


# ! hacer este para un proceso multidimensional, para unidimensional no tiene sentido
# Case of study:
#  - OneDimensional -> D == 1
#  - NonDiagonalNoise with M == D but NonDiagonal matrix
#  - Out Of Place -> IIP = false
#  - Without Interest Rates
# SDEs:
#   dx(t) = r_x ⋅ x(t) ⋅ dt + σ_x ⋅ x(t) ⋅ dW_x(t)



# Case of study:
#  - MultiDimensional -> D > 1
#  - DiagonalNoise (all subsystems have DiagonalNoise) -> M = D
#  - Out Of Place -> IIP = false
#  - Without Interest Rates
# SDEs:
#   dx(t) = r_x ⋅ x(t) ⋅ dt + σ_x ⋅ x(t) ⋅ dW_x(t)
#   dY(t) = r_Y ⋅ Y(t) ⋅ dt + σ_Y ⋅ Y(t) ⋅ dW_Y(t)
#   dz(t) = r_z ⋅ z(t) ⋅ dt + σ_z ⋅ z(t) ⋅ dW_z(t)
#   dK(t) = r_K ⋅ K(t) ⋅ dt + σ_K ⋅ K(t) ⋅ dW_K(t)
# Note that each diffusion must yield to a vector since we have DiagonalNoise.

# The first step is to describe each dynamics present in the system using its initial
# condition, noise type, correlations and initial time.

# x is a one dimensional security process with DiagonalNoise
x0 = 0.1 # SVector(0.1)
x_dynamics = SystemDynamics(x0; t0=zero(x0), m=DiagonalNoise(1))

# Y is a multi dimensional security process with DiagonalNoise
NY = 3
Y0 = @SVector [0.2, 0.3, 0.4]
Y_dynamics = SystemDynamics(Y0, t0=zero(eltype(Y0)), m=DiagonalNoise(NY))

# z is a one dimensional security process with DiagonalNoise (as x)
z0 = 0.5 # SVector(0.5)
z_dynamics = SystemDynamics(z0; t0=zero(x0), m=DiagonalNoise(1))

# K is a multi dimensional security process with DiagonalNoise (as Y)
NK = 2
K0 = @SVector [0.6, 0.7]
K_dynamics = SystemDynamics(K0, t0=zero(eltype(K0)), m=DiagonalNoise(NK))

# Once we have all the dynamics that compose the system, we can wrap them in a container
# and compute some relevant system attributes
dynamics = (
    x_dynamics = x_dynamics,
    Y_dynamics = Y_dynamics,
    z_dynamics = z_dynamics,
    K_dynamics = K_dynamics
)

r_x = 0.01
r_Y = @SMatrix[0.01 0.02 0.03; 0.04 0.05 0.06; 0.07 0.08 0.09]
r_z = 0.02
r_K = @SMatrix[0.01 0.02; 0.03 0.04]

σ_x = r_x
σ_Y = r_Y
σ_z = r_z
σ_K = r_K

# Podria haber definido todo lo anterior dentro de params
params = (
    r_x = r_x,
    r_Y = r_Y,
    r_z = r_z,
    r_K = r_K,
    σ_x = σ_x,
    σ_Y = σ_Y,
    σ_z = σ_z,
    σ_K = σ_K,
    dynamics...
)

function system_μ(u, p, t)
    @unpack r_x, r_Y, r_x, r_K, σ_x, σ_Y, σ_x, σ_K, x_dynamics, Y_dynamics, z_dynamics, K_dynamics = p

    x = Security{false,1,1,true}(u, 1, t)
    Y = Security{false,3,3,true}(u, 2:4, t)
    z = Security{false,1,1,true}(u, 5, t)
    K = Security{false,2,2,true}(u, 6:7, t)

    dx = r_x * x(t)
    dY = r_Y * Y(t)
    dz = r_z * z(t)
    dK = r_K * K(t)

    # return vcat(dx, dY, dz, dK) # aun no funciona bien, ya que traduce a Array
    return vcat(SVector(dx), dY, SVector(dz), dK)
end

function system_σ(u, p, t)
    @unpack r_x, r_Y, r_x, r_K, σ_x, σ_Y, σ_x, σ_K, x_dynamics, Y_dynamics, z_dynamics, K_dynamics = p

    x = Security{false,1,1,true}(u, 1, t)
    Y = Security{false,3,3,true}(u, 2:4, t)
    z = Security{false,1,1,true}(u, 5, t)
    K = Security{false,2,2,true}(u, 6:7, t)

    dx = σ_x * x(t)
    dY = σ_Y * Y(t)
    dz = σ_z * z(t)
    dK = σ_K * K(t)

    # return vcat(dx, dY, dz, dK) # aun no funciona bien, ya que traduce a Array
    return vcat(SVector(dx), dY, SVector(dz), dK)
end

ds = DynamicalSystem(system_μ, system_σ, params)

@test isinplace(ds) == false
@test dimension(ds) == 7
@test noise_dimension(ds) == 7
@test diagonalnoise(ds) == true
#! agregar muchos mas tests

# Case of study:
#  - MultiDimensional -> D > 1
#  - All subsystems have ScalarNoise -> System ends up with NonDiagonalNoise
#  - Out Of Place -> IIP = false
#  - Without Interest Rates
# SDEs:
#   dx(t) = r_x ⋅ x(t) ⋅ dt + σ_x ⋅ x(t) ⋅ dW₁(t)
#   dY(t) = r_Y ⋅ Y(t) ⋅ dt + σ_Y ⋅ Y(t) ⋅ dW₂(t)
#   dz(t) = r_z ⋅ z(t) ⋅ dt + σ_z ⋅ z(t) ⋅ dW₃(t)
#   dK(t) = r_K ⋅ K(t) ⋅ dt + σ_K ⋅ K(t) ⋅ dW₄(t)
# Note that each diffusion must yield to a scalar since we have ScalarNoise.

# The first step is to describe each dynamics present in the system using its initial
# condition, noise type, correlations and initial time.

# x is a one dimensional security process with DiagonalNoise
x0 = 0.1 # SVector(0.1)
x_dynamics = SystemDynamics(x0; t0=zero(x0), m=ScalarNoise())

# Y is a multi dimensional security process with DiagonalNoise
NY = 3
Y0 = @SVector [0.2, 0.3, 0.4]
Y_dynamics = SystemDynamics(Y0, t0=zero(eltype(Y0)), m=ScalarNoise())

# z is a one dimensional security process with DiagonalNoise (as x)
z0 = 0.5 # SVector(0.5)
z_dynamics = SystemDynamics(z0; t0=zero(x0), m=ScalarNoise())

# K is a multi dimensional security process with DiagonalNoise (as Y)
NK = 2
K0 = @SVector [0.6, 0.7]
K_dynamics = SystemDynamics(K0, t0=zero(eltype(K0)), m=ScalarNoise())

# Once we have all the dynamics that compose the system, we can wrap them in a container
# and compute some relevant system attributes
dynamics = (
    x_dynamics = x_dynamics,
    Y_dynamics = Y_dynamics,
    z_dynamics = z_dynamics,
    K_dynamics = K_dynamics
)

sys_attrs = DynamicalSystemAttributes(dynamics)

r_x = 0.01
r_Y = @SMatrix[0.01 0.02 0.03; 0.04 0.05 0.06; 0.07 0.08 0.09]
r_z = 0.02
r_K = @SMatrix[0.01 0.02; 0.03 0.04]

σ_x = r_x
σ_Y = r_Y
σ_z = r_z
σ_K = r_K

# Podria haber definido todo lo anterior dentro de params
params = (
    r_x = r_x,
    r_Y = r_Y,
    r_z = r_z,
    r_K = r_K,
    σ_x = σ_x,
    σ_Y = σ_Y,
    σ_z = σ_z,
    σ_K = σ_K,
    dynamics...
)

function system_μ(u, p, t)
    @unpack r_x, r_Y, r_x, r_K, σ_x, σ_Y, σ_x, σ_K, x_dynamics, Y_dynamics, z_dynamics, K_dynamics = p

    x = Security{false,1,1,true}(u, 1, t)
    Y = Security{false,3,3,true}(u, 2:4, t)
    z = Security{false,1,1,true}(u, 5, t)
    K = Security{false,2,2,true}(u, 6:7, t)

    dx = r_x * x(t)
    dY = r_Y * Y(t)
    dz = r_z * z(t)
    dK = r_K * K(t)

    # return vcat(dx, dY, dz, dK) # aun no funciona bien, ya que traduce a Array
    return vcat(SVector(dx), dY, SVector(dz), dK)
end

function system_σ(u, p, t)
    @unpack r_x, r_Y, r_x, r_K, σ_x, σ_Y, σ_x, σ_K, x_dynamics, Y_dynamics, z_dynamics, K_dynamics = p

    x = Security{false,1,1,true}(u, 1, t)
    Y = Security{false,3,3,true}(u, 2:4, t)
    z = Security{false,1,1,true}(u, 5, t)
    K = Security{false,2,2,true}(u, 6:7, t)

    dx = σ_x * x(t)
    dY = σ_Y * Y(t)
    dz = σ_z * z(t)
    dK = σ_K * K(t)

    # return vcat(dx, dY, dz, dK) # aun no funciona bien, ya que traduce a Array
    return vcat(SVector(dx), dY, SVector(dz), dK)
end

ds = DynamicalSystem(system_μ, system_σ, params)

@test isinplace(ds) == false
@test dimension(ds) == 7
@test noise_dimension(ds) == 7
@test diagonalnoise(ds) == true