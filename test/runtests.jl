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
        @test isa(attrs.noise_rate_prototype, SMatrix{(D,M)...,T})
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
        @test attrs.noise == WienerProcess(attrs.t0, @SVector(zeros(M)), @SVector(zeros(M)))
        @test isnothing(attrs.noise_rate_prototype) # diffeq will use u0 as prototype!

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        @test isinplace(ds) == IIP
        @test dimension(ds) == D
        @test noise_dimension(ds) == M
        @test diagonalnoise(ds) == DN
        @test UniversalDynamics.eltype(ds) == T
    end

    @testset "MultiDimensional + DiagonalNoise" begin
        IIP = false
        D = 3
        M = 3
        DN = true
        T = Float64

        t0 = zero(T)
        x0 = @SVector rand(D)
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
        @test state(attrs) == x0
        @test cor(attrs) == I
        @test isnothing(attrs.noise)
        @test isnothing(attrs.noise_rate_prototype)

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        @test isinplace(ds) == IIP
        @test dimension(ds) == D
        @test noise_dimension(ds) == M
        @test diagonalnoise(ds) == DN
        @test UniversalDynamics.eltype(ds) == T
    end

    @testset "MultiDimensional + NonDiagonalNoise" begin
        IIP = false
        D = 3
        M = 5
        DN = false
        T = Float64

        t0 = zero(T)
        x0 = @SVector rand(D)
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
        @test isa(attrs.noise_rate_prototype, SMatrix{(D,M)...,T})
        @test attrs.noise_rate_prototype == @SMatrix(ones(D, M))

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, )) # dummy `f` and `g`
        @test isinplace(ds) == IIP
        @test dimension(ds) == D
        @test noise_dimension(ds) == M
        @test diagonalnoise(ds) == DN
        @test UniversalDynamics.eltype(ds) == T
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

        x0 = 0.1
        x_dynamics = SystemDynamics(x0; t0=t0, noise=ScalarNoise())

        @test isinplace(x_dynamics) == IIP
        @test dimension(x_dynamics) == Dx
        @test noise_dimension(x_dynamics) == Mx
        @test diagonalnoise(x_dynamics) == DNx
        @test UniversalDynamics.eltype(x_dynamics) == T
        @test initialtime(x_dynamics) == t0
        @test state(x_dynamics) == x0
        @test cor(x_dynamics) == Diagonal(SVector{Mx,T}(ones(Mx))) == I

        y0 = @SVector rand(Dy)
        y_dynamics = SystemDynamics(y0; t0=t0, noise=NonDiagonalNoise(My))

        @test isinplace(y_dynamics) == IIP
        @test dimension(y_dynamics) == Dy
        @test noise_dimension(y_dynamics) == My
        @test diagonalnoise(y_dynamics) == DNy
        @test UniversalDynamics.eltype(y_dynamics) == T
        @test initialtime(y_dynamics) == t0
        @test state(y_dynamics) == y0
        @test cor(y_dynamics) == Diagonal(SVector{My,T}(ones(My))) == I
        @test UniversalDynamics.gprototype(y_dynamics) == @SVector(ones(Dy))

        attrs = DynamicalSystemAttributes((x_dynamics, y_dynamics, ))
        @test initialtime(attrs) == t0
        @test state(attrs) == vcat(SVector(x0), y0)
        @test cor(attrs) == I
        @test attrs.noise == WienerProcess(attrs.t0, @SVector(zeros(M)), @SVector(zeros(M)))
        @test isa(attrs.noise_rate_prototype, SMatrix{(D,M)...,T})
        @test attrs.noise_rate_prototype == @SMatrix [1. 0.; 0. 1.; 0. 1.; 0. 1.]

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, y_dynamics, )) # dummy `f` and `g`
        @test isinplace(ds) == IIP
        @test dimension(ds) == D
        @test noise_dimension(ds) == M
        @test diagonalnoise(ds) == DN
        @test UniversalDynamics.eltype(ds) == T
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

        x0 = 0.1
        x_dynamics = SystemDynamics(x0; t0=t0, noise=DiagonalNoise(Mx))

        @test isinplace(x_dynamics) == IIP
        @test dimension(x_dynamics) == Dx
        @test noise_dimension(x_dynamics) == Mx
        @test diagonalnoise(x_dynamics) == DNx
        @test UniversalDynamics.eltype(x_dynamics) == T
        @test initialtime(x_dynamics) == t0
        @test state(x_dynamics) == x0
        @test cor(x_dynamics) == Diagonal(SVector{Mx,T}(ones(Mx))) == I

        y0 = @SVector rand(Dy)
        y_dynamics = SystemDynamics(y0; t0=t0, noise=DiagonalNoise(My))

        @test isinplace(y_dynamics) == IIP
        @test dimension(y_dynamics) == Dy
        @test noise_dimension(y_dynamics) == My
        @test diagonalnoise(y_dynamics) == DNy
        @test UniversalDynamics.eltype(y_dynamics) == T
        @test initialtime(y_dynamics) == t0
        @test state(y_dynamics) == y0
        @test cor(y_dynamics) == Diagonal(SVector{My,T}(ones(My))) == I
        @test UniversalDynamics.gprototype(y_dynamics) == @SVector(ones(Dy))

        attrs = DynamicalSystemAttributes((x_dynamics, y_dynamics, ))
        @test initialtime(attrs) == t0
        @test state(attrs) == vcat(SVector(x0), y0)
        @test cor(attrs) == I
        @test isnothing(attrs.noise)
        @test isnothing(attrs.noise_rate_prototype)

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, y_dynamics, )) # dummy `f` and `g`
        @test isinplace(ds) == IIP
        @test dimension(ds) == D
        @test noise_dimension(ds) == M
        @test diagonalnoise(ds) == DN
        @test UniversalDynamics.eltype(ds) == T
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

        @test isinplace(x_dynamics) == IIP
        @test dimension(x_dynamics) == Dx
        @test noise_dimension(x_dynamics) == Mx
        @test diagonalnoise(x_dynamics) == DNx
        @test UniversalDynamics.eltype(x_dynamics) == T
        @test initialtime(x_dynamics) == t0
        @test state(x_dynamics) == SVector(x0)
        @test cor(x_dynamics) == Diagonal(SVector{Mx,T}(ones(Mx))) == I
        @test UniversalDynamics.gprototype(x_dynamics) == @SMatrix(ones(Dx, Mx))

        y0 = @SVector rand(Dy)
        y_dynamics = SystemDynamics(y0; t0=t0, noise=NonDiagonalNoise(My))

        @test isinplace(y_dynamics) == IIP
        @test dimension(y_dynamics) == Dy
        @test noise_dimension(y_dynamics) == My
        @test diagonalnoise(y_dynamics) == DNy
        @test UniversalDynamics.eltype(y_dynamics) == T
        @test initialtime(y_dynamics) == t0
        @test state(y_dynamics) == y0
        @test cor(y_dynamics) == Diagonal(SVector{My,T}(ones(My))) == I
        @test UniversalDynamics.gprototype(y_dynamics) == @SMatrix(ones(Dy, My))

        attrs = DynamicalSystemAttributes((x_dynamics, y_dynamics, ))
        @test initialtime(attrs) == t0
        @test state(attrs) == vcat(x0, y0)
        @test cor(attrs) == I
        @test attrs.noise == WienerProcess(attrs.t0, @SVector(zeros(M)), @SVector(zeros(M)))
        @test attrs.noise_rate_prototype == @SMatrix [1. 1. 1. 1. 0. 0. 0. 0. 0.; 0. 0. 0. 0. 1. 1. 1. 1. 1.; 0. 0. 0. 0. 1. 1. 1. 1. 1.; 0. 0. 0. 0. 1. 1. 1. 1. 1.]

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, y_dynamics, )) # dummy `f` and `g`
        @test isinplace(ds) == IIP
        @test dimension(ds) == D
        @test noise_dimension(ds) == M
        @test diagonalnoise(ds) == DN
        @test UniversalDynamics.eltype(ds) == T
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

        x0 = 0.1
        x_dynamics = SystemDynamics(x0; t0=t0, noise=ScalarNoise())

        @test isinplace(x_dynamics) == IIP
        @test dimension(x_dynamics) == Dx
        @test noise_dimension(x_dynamics) == Mx
        @test diagonalnoise(x_dynamics) == DNx
        @test UniversalDynamics.eltype(x_dynamics) == T
        @test initialtime(x_dynamics) == t0
        @test state(x_dynamics) == x0
        @test cor(x_dynamics) == Diagonal(SVector{Mx,T}(ones(Mx))) == I

        y0 = @SVector rand(Dy)
        y_dynamics = SystemDynamics(y0; t0=t0, noise=DiagonalNoise(My))

        @test isinplace(y_dynamics) == IIP
        @test dimension(y_dynamics) == Dy
        @test noise_dimension(y_dynamics) == My
        @test diagonalnoise(y_dynamics) == DNy
        @test UniversalDynamics.eltype(y_dynamics) == T
        @test initialtime(y_dynamics) == t0
        @test state(y_dynamics) == y0
        @test cor(y_dynamics) == Diagonal(SVector{My,T}(ones(My))) == I
        @test UniversalDynamics.gprototype(y_dynamics) == @SVector(ones(Dy))

        z0 = @SVector rand(Dz)
        z_dynamics = SystemDynamics(z0; t0=t0, noise=NonDiagonalNoise(Mz))

        @test isinplace(z_dynamics) == IIP
        @test dimension(z_dynamics) == Dz
        @test noise_dimension(z_dynamics) == Mz
        @test diagonalnoise(z_dynamics) == DNz
        @test UniversalDynamics.eltype(z_dynamics) == T
        @test initialtime(z_dynamics) == t0
        @test state(z_dynamics) == z0
        @test cor(z_dynamics) == Diagonal(SVector{Mz,T}(ones(Mz))) == I
        @test UniversalDynamics.gprototype(z_dynamics) == @SMatrix(ones(Dz, Mz))

        attrs = DynamicalSystemAttributes((x_dynamics, y_dynamics, z_dynamics))
        @test initialtime(attrs) == t0
        @test state(attrs) == vcat(SVector(x0), y0, z0)
        @test cor(attrs) == I
        @test attrs.noise == WienerProcess(attrs.t0, @SVector(zeros(M)), @SVector(zeros(M)))
        @test isa(attrs.noise_rate_prototype, SMatrix{(D,M)...,T})
        @test attrs.noise_rate_prototype == @SMatrix [1. 0. 0. 0. 0.; 0. 1. 0. 0. 0.; 0. 0. 1. 1. 1.]

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, y_dynamics, z_dynamics)) # dummy `f` and `g`
        @test isinplace(ds) == IIP
        @test dimension(ds) == D
        @test noise_dimension(ds) == M
        @test diagonalnoise(ds) == DN
        @test UniversalDynamics.eltype(ds) == T
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

        @test isinplace(x_dynamics) == IIP
        @test dimension(x_dynamics) == Dx
        @test noise_dimension(x_dynamics) == Mx
        @test diagonalnoise(x_dynamics) == DNx
        @test UniversalDynamics.eltype(x_dynamics) == T
        @test initialtime(x_dynamics) == t0
        @test state(x_dynamics) == x0
        @test cor(x_dynamics) == Diagonal(SVector{Mx,T}(ones(Mx))) == I

        y0 = @SVector rand(Dy)
        y_dynamics = SystemDynamics(y0; t0=t0, noise=DiagonalNoise(My))

        @test isinplace(y_dynamics) == IIP
        @test dimension(y_dynamics) == Dy
        @test noise_dimension(y_dynamics) == My
        @test diagonalnoise(y_dynamics) == DNy
        @test UniversalDynamics.eltype(y_dynamics) == T
        @test initialtime(y_dynamics) == t0
        @test state(y_dynamics) == y0
        @test cor(y_dynamics) == Diagonal(SVector{My,T}(ones(My))) == I
        @test UniversalDynamics.gprototype(y_dynamics) == @SVector(ones(Dy))

        z0 = @SVector rand(Dz)
        z_dynamics = SystemDynamics(z0; t0=t0, noise=NonDiagonalNoise(Mz))

        @test isinplace(z_dynamics) == IIP
        @test dimension(z_dynamics) == Dz
        @test noise_dimension(z_dynamics) == Mz
        @test diagonalnoise(z_dynamics) == DNz
        @test UniversalDynamics.eltype(z_dynamics) == T
        @test initialtime(z_dynamics) == t0
        @test state(z_dynamics) == z0
        @test cor(z_dynamics) == Diagonal(SVector{Mz,T}(ones(Mz))) == I
        @test UniversalDynamics.gprototype(z_dynamics) == @SMatrix(ones(Dz, Mz))

        attrs = DynamicalSystemAttributes((x_dynamics, y_dynamics, z_dynamics))
        @test initialtime(attrs) == t0
        @test state(attrs) == vcat(x0, y0, z0)
        @test cor(attrs) == I
        @test attrs.noise == WienerProcess(attrs.t0, @SVector(zeros(M)), @SVector(zeros(M)))
        @test isa(attrs.noise_rate_prototype, SMatrix{(D,M)...,T})
        @test attrs.noise_rate_prototype == @SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 1.0 1.0 1.0; 0.0 0.0 0.0 0.0 1.0 1.0 1.0; 0.0 0.0 0.0 0.0 1.0 1.0 1.0; 0.0 0.0 0.0 0.0 1.0 1.0 1.0]

        ds = DynamicalSystem(nothing, nothing, (x_dynamics, y_dynamics, z_dynamics)) # dummy `f` and `g`
        @test isinplace(ds) == IIP
        @test dimension(ds) == D
        @test noise_dimension(ds) == M
        @test diagonalnoise(ds) == DN
        @test UniversalDynamics.eltype(ds) == T
    end
end