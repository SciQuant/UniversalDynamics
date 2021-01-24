using UniversalDynamics
using Test

using DiffEqNoiseProcess
using LinearAlgebra
using StaticArrays


function test_dynamics(d, IIP, D, M, DN, T, t0, x0, ρ, noise, noise_rate_prototype)
    @test isinplace(d) == IIP
    @test dimension(d) == D
    @test noise_dimension(d) == M
    @test diagonalnoise(d) == DN
    @test UniversalDynamics.eltype(d) == T

    @test initialtime(d) == t0 && typeof(initialtime(d)) == typeof(t0)
    @test state(d) == x0 && typeof(state(d)) == typeof(x0)
    @test cor(d) == ρ #&& typeof(cor(d)) == ρ
    @test UniversalDynamics.noise(d) == noise
    @test typeof(UniversalDynamics.noise(d)) == typeof(noise)
    @test UniversalDynamics.noise_rate_prototype(d) == noise_rate_prototype
    @test typeof(UniversalDynamics.noise_rate_prototype(d)) == typeof(noise_rate_prototype)

    return nothing
end

@testset "Out of Place" begin

    @testset "OneDimensional + ScalarNoise" begin
        IIP, D, M, DN, T = false, 1, 1, true, Float64

        t0 = zero(T)

        x0 = T(0.1) # SVector(0.1)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=ScalarNoise())
        test_dynamics(
            x_dynamics, IIP, D, M, DN, T,
            t0, x0, I, nothing, @SVector(ones(D))
        )

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        test_dynamics(
            ds,
            IIP, D, M, DN, T,
            t0, SVector(x0), I, nothing, nothing
        )
    end

    # should yield to same results as the "OneDimensional + ScalarNoise" test
    @testset "OneDimensional + DiagonalNoise" begin
        IIP, D, M, DN, T = false, 1, 1, true, Float64

        t0 = zero(T)

        x0 = T(0.1) # SVector(0.1)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=DiagonalNoise(M))
        test_dynamics(x_dynamics, IIP, D, M, DN, T, t0, x0, I, nothing, @SVector(ones(D)))

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        test_dynamics(ds, IIP, D, M, DN, T, t0, SVector(x0), I, nothing, nothing)
    end

    @testset "OneDimensional + NonDiagonalNoise" begin
        IIP, D, M, DN, T = false, 1, 3, false, Float64

        t0 = zero(T)

        x0 = T(0.1)
        @test_throws ArgumentError SystemDynamics(x0; t0=t0, noise=NonDiagonalNoise(M))
        x0 = SVector{1,T}(0.1) # must be a Vector
        x_dynamics = SystemDynamics(x0; t0=t0, noise=NonDiagonalNoise(M))
        test_dynamics(
            x_dynamics,
            IIP, D, M, DN, T,
            t0, x0, I,
            WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
            @SMatrix(ones(D, M))
        )

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        test_dynamics(
            ds,
            IIP, D, M, DN, T,
            t0, x0, I,
            WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
            @SMatrix(ones(D, M))
        )
    end

    @testset "MultiDimensional + ScalarNoise" begin
        IIP, D, M, DN, T = false, 3, 1, false, Float64

        t0 = zero(T)

        x0 = @SVector rand(D)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=ScalarNoise())
        test_dynamics(
            x_dynamics,
            IIP, D, M, DN, T,
            t0, x0, I,
            WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
            @SVector(ones(D))
        )

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        test_dynamics(
            ds,
            IIP, D, M, DN, T,
            t0, x0, I,
            WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
            nothing # diffeq will use u0 as gprototype!
        )
    end

    @testset "MultiDimensional + DiagonalNoise" begin
        IIP, D, M, DN, T = false, 3, 3, true, Float64

        t0 = zero(T)

        x0 = @SVector rand(D)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=DiagonalNoise(M))
        test_dynamics(
            x_dynamics,
            IIP, D, M, DN, T,
            t0, x0, I, nothing, @SVector(ones(D))
        )

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        test_dynamics(
            ds,
            IIP, D, M, DN, T,
            t0, x0, I, nothing, nothing
        )
    end

    @testset "MultiDimensional + NonDiagonalNoise" begin
        IIP, D, M, DN, T = false, 3, 5, false, Float64

        t0 = zero(T)

        x0 = @SVector rand(D)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=NonDiagonalNoise(M))
        test_dynamics(
            x_dynamics,
            IIP, D, M, DN, T,
            t0, x0, I,
            WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
            @SMatrix(ones(D, M))
        )

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        test_dynamics(
            ds,
            IIP, D, M, DN, T,
            t0, x0, I,
            WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
            @SMatrix(ones(D, M))
        )
    end

    @testset "Mixed One and Multi Dimensional + ScalarNoises" begin
        IIP = false

        Dx = 1
        Mx = 1
        DNx = true

        Dy = 3
        My = 1
        DNy = false

        D = Dx + Dy
        M = Mx + My
        DN = false

        T = Float64

        t0 = zero(T)

        x0 = T(0.1)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=ScalarNoise())
        test_dynamics(
            x_dynamics,
            IIP, Dx, Mx, DNx, T,
            t0, x0, I, nothing, @SVector(ones(Dx))
        )

        y0 = @SVector rand(Dy)
        y_dynamics = SystemDynamics(y0; t0=t0, noise=ScalarNoise())
        test_dynamics(
            y_dynamics,
            IIP, Dy, My, DNy, T,
            t0, y0, I,
            WienerProcess(t0, @SVector(zeros(My)), @SVector(zeros(My))),
            @SVector(ones(Dy))
        )

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, y_dynamics)) # dummy `f` and `g`
        test_dynamics(
            ds,
            IIP, D, M, DN, T,
            t0, vcat(SVector(x0), y0), I,
            WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
            @SMatrix [1. 0.; 0. 1.; 0. 1.; 0. 1.]
        )
    end

    @testset "Mixed One and Multi Dimensional + DiagonalNoises" begin
        IIP = false

        Dx = 1
        Mx = 1
        DNx = true

        Dy = 3
        My = 3
        DNy = true

        D = Dx + Dy
        M = Mx + My
        DN = true

        T = Float64

        t0 = zero(T)

        x0 = T(0.1)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=DiagonalNoise(Mx))
        test_dynamics(
            x_dynamics,
            IIP, Dx, Mx, DNx, T,
            t0, x0, I, nothing, @SVector(ones(Dx))
        )

        y0 = @SVector rand(Dy)
        y_dynamics = SystemDynamics(y0; t0=t0, noise=DiagonalNoise(My))
        test_dynamics(
            y_dynamics,
            IIP, Dy, My, DNy, T,
            t0, y0, I, nothing, @SVector(ones(Dy))
        )

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, y_dynamics)) # dummy `f` and `g`
        test_dynamics(
            ds,
            IIP, D, M, DN, T,
            t0, vcat(SVector(x0), y0), I, nothing, nothing
        )
    end

    @testset "Mixed One and Multi Dimensional + NonDiagonalNoises" begin
        IIP = false

        Dx = 1
        Mx = 4
        DNx = false

        Dy = 3
        My = 5
        DNy = false

        D = Dx + Dy
        M = Mx + My
        DN = false

        T = Float64

        t0 = zero(T)

        x0 = SVector(0.1)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=NonDiagonalNoise(Mx))
        test_dynamics(
            x_dynamics,
            IIP, Dx, Mx, DNx, T,
            t0, x0, I,
            WienerProcess(t0, @SVector(zeros(Mx)), @SVector(zeros(Mx))),
            @SMatrix(ones(Dx, Mx))
        )

        y0 = @SVector rand(Dy)
        y_dynamics = SystemDynamics(y0; t0=t0, noise=NonDiagonalNoise(My))
        test_dynamics(
            y_dynamics,
            IIP, Dy, My, DNy, T,
            t0, y0, I,
            WienerProcess(t0, @SVector(zeros(My)), @SVector(zeros(My))),
            @SMatrix(ones(Dy, My))
        )

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, y_dynamics)) # dummy `f` and `g`
        test_dynamics(
            ds,
            IIP, D, M, DN, T,
            t0, vcat(SVector(x0), y0), I,
            WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
            @SMatrix [1. 1. 1. 1. 0. 0. 0. 0. 0.; 0. 0. 0. 0. 1. 1. 1. 1. 1.; 0. 0. 0. 0. 1. 1. 1. 1. 1.; 0. 0. 0. 0. 1. 1. 1. 1. 1.]
        )
    end

    @testset "OneDimensional + Mixed Noises" begin
        IIP = false

        # ScalarNoise
        Dx = 1
        Mx = 1
        DNx = true

        # DiagonalNoise (for the one dimensional case it the same as ScalarNoise)
        Dy = 1
        My = 1
        DNy = true

        # NonDiagonalNoise
        Dz = 1
        Mz = 3
        DNz = false

        D = Dx + Dy + Dz
        M = Mx + My + Mz
        DN = false

        T = Float64

        t0 = zero(T)

        x0 = T(0.1)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=ScalarNoise())
        test_dynamics(
            x_dynamics,
            IIP, Dx, Mx, DNx, T,
            t0, x0, I,
            nothing,
            @SVector(ones(Dx))
        )

        y0 = @SVector rand(Dy)
        y_dynamics = SystemDynamics(y0; t0=t0, noise=DiagonalNoise(My))
        test_dynamics(
            y_dynamics,
            IIP, Dy, My, DNy, T,
            t0, y0, I,
            nothing,
            @SVector(ones(Dy))
        )

        z0 = @SVector rand(Dz)
        z_dynamics = SystemDynamics(z0; t0=t0, noise=NonDiagonalNoise(Mz))
        test_dynamics(
            z_dynamics,
            IIP, Dz, Mz, DNz, T,
            t0, z0, I,
            WienerProcess(t0, @SVector(zeros(Mz)), @SVector(zeros(Mz))),
            @SMatrix(ones(Dz, Mz))
        )

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, y_dynamics, z_dynamics)) # dummy `f` and `g`
        test_dynamics(
            ds,
            IIP, D, M, DN, T,
            t0, vcat(SVector(x0), y0, z0), I,
            WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
            @SMatrix [1. 0. 0. 0. 0.; 0. 1. 0. 0. 0.; 0. 0. 1. 1. 1.]
        )
    end

    @testset "MultiDimensional + Mixed Noises" begin
        IIP = false

        # ScalarNoise
        Dx = 2
        Mx = 1
        DNx = false

        # DiagonalNoise (for the one dimensional case it the same as ScalarNoise)
        Dy = 3
        My = 3
        DNy = true

        # NonDiagonalNoise
        Dz = 4
        Mz = 3
        DNz = false

        D = Dx + Dy + Dz
        M = Mx + My + Mz
        DN = false

        T = Float64

        t0 = zero(T)

        x0 = @SVector rand(Dx)
        x_dynamics = SystemDynamics(x0; t0=t0, noise=ScalarNoise())
        test_dynamics(
            x_dynamics,
            IIP, Dx, Mx, DNx, T,
            t0, x0, I,
            WienerProcess(t0, @SVector(zeros(Mx)), @SVector(zeros(Mx))),
            @SVector(ones(Dx))
        )

        y0 = @SVector rand(Dy)
        y_dynamics = SystemDynamics(y0; t0=t0, noise=DiagonalNoise(My))
        test_dynamics(
            y_dynamics,
            IIP, Dy, My, DNy, T,
            t0, y0, I,
            nothing,
            @SVector(ones(Dy))
        )

        z0 = @SVector rand(Dz)
        z_dynamics = SystemDynamics(z0; t0=t0, noise=NonDiagonalNoise(Mz))
        test_dynamics(
            z_dynamics,
            IIP, Dz, Mz, DNz, T,
            t0, z0, I,
            WienerProcess(t0, @SVector(zeros(Mz)), @SVector(zeros(Mz))),
            @SMatrix(ones(Dz, Mz))
        )

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, y_dynamics, z_dynamics)) # dummy `f` and `g`
        test_dynamics(
            ds,
            IIP, D, M, DN, T,
            t0, vcat(SVector(x0), y0, z0), I,
            WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
            @SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 1.0 1.0; 0.0 0.0 0.0 0.0 1.0 1.0 1.0; 0.0 0.0 0.0 0.0 1.0 1.0 1.0; 0.0 0.0 0.0 0.0 1.0 1.0 1.0]
        )
    end
end

r0 = ones(1)
κ = θ = Σ = α = β = one
r = OneFactorAffineModelDynamics(r0, κ, θ, Σ, α, β)

x0 = rand(3)
κ = θ = α = β = Σ = ξ₀ = ξ₁ = one
Σd(t) = Diagonal(rand(3))
x = MultiFactorAffineModelDynamics(x0, κ, θ, Σd, α, β, ξ₀, ξ₁)

Dy = 3
My = 4
y0 = rand(Dy)
y = SystemDynamics(y0; noise=NonDiagonalNoise(My))

dynamics = (r, x, y)
ds = DynamicalSystem(dynamics)


include("../../UniversalMonteCarlo/test/DaiSingletonParameters_A3_1.jl")
N = 3
# x0 = @SVector [υ₀, θ₀, r₀]
x0 = @SVector ones(N)

ξ₀(t) = zero(t) # ξ₀ = zero
ξ₁(t) = @SVector [0, 0, 1]

ϰ(t) = @SMatrix([
    μ     0 0
    0     ν 0
    κ_rυ -κ κ
])
θ(t) = @SVector [ῡ, θ̄, θ̄ ]
Σ(t) = @SMatrix [
    η           0    0
    η * σ_θυ    1 σ_θr
    η * σ_rυ σ_rθ    1
]

α(t) = @SVector [0, ζ^2, α_r]
β(t) = @SMatrix [
    1 0 0
    β_θ 0 0
    1 0 0
]

x = MultiFactorAffineModelDynamics(x0, ϰ, θ, Σ, α, β, ξ₀, ξ₁)
B = SystemDynamics(0.)

ds = DynamicalSystem((x, B))

IR = FixedIncomeSecurities(x, ds.securities[1], ds.securities[2])

IR.P(0.15, 0.15)

function drift(u, p, t)
    @unpack security_x = p

    x = remake(security_x, u)

    # ambos funcionan
    x()
    x(t)
end

function drift!(du, u, p, t)
    @unpack security_x = p

    x = remake(security_x, u, du)

    x()
    x(t)
    x.dx # ver de agregar σ(x) or μ(x)
end