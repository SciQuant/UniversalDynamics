
function test_dynamics(d, IIP, D, M, DN, T, t0, x0, ρ, noise, noise_rate_prototype)
    @test isinplace(d) == IIP
    @test dimension(d) == D
    @test noise_dimension(d) == M
    @test diagonalnoise(d) == DN
    @test eltype(d) == T

    @test get_t0(d) == t0 && typeof(get_t0(d)) == typeof(t0)
    @test get_state(d) == x0 && typeof(get_state(d)) == typeof(x0)
    @test get_cor(d) == ρ #&& typeof(get_cor(d)) == ρ
    @test get_noise(d) == noise
    @test typeof(get_noise(d)) == typeof(noise)
    @test get_noise_rate_prototype(d) == noise_rate_prototype
    @test typeof(get_noise_rate_prototype(d)) == typeof(noise_rate_prototype)

    return nothing
end

@testset "OneDimensional + ScalarNoise" begin
    IIP, D, M, DN, T = false, 1, 1, true, Float64

    t0 = zero(T)

    x0 = T(0.1) # SVector(0.1)
    x_dynamics = SystemDynamics(x0; t0=t0, noise=ScalarNoise())
    test_dynamics(
        x_dynamics, IIP, D, M, DN, T,
        t0, x0, I, nothing, @SVector(ones(D))
    )

    ds = DynamicalSystem([:x => x_dynamics])
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

    ds = DynamicalSystem([:x => x_dynamics])
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

    ds = DynamicalSystem([:x => x_dynamics])
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

    ds = DynamicalSystem([:x => x_dynamics])
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

    ds = DynamicalSystem([:x => x_dynamics])
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

    ds = DynamicalSystem([:x => x_dynamics])
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

    ds = DynamicalSystem([:x => x_dynamics, :y => y_dynamics])
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

    ds = DynamicalSystem([:x => x_dynamics, :y => y_dynamics])
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

    ds = DynamicalSystem([:x => x_dynamics, :y => y_dynamics])
    test_dynamics(
        ds,
        IIP, D, M, DN, T,
        t0, vcat(SVector(x0), y0), I,
        WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
        @SMatrix [1. 1. 1. 1. 0. 0. 0. 0. 0.
                  0. 0. 0. 0. 1. 1. 1. 1. 1.
                  0. 0. 0. 0. 1. 1. 1. 1. 1.
                  0. 0. 0. 0. 1. 1. 1. 1. 1.]
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

    ds = DynamicalSystem([:x => x_dynamics, :y => y_dynamics, :z => z_dynamics])
    test_dynamics(
        ds,
        IIP, D, M, DN, T,
        t0, vcat(SVector(x0), y0, z0), I,
        WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
        @SMatrix [1. 0. 0. 0. 0.
                  0. 1. 0. 0. 0.
                  0. 0. 1. 1. 1.]
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

    ds = DynamicalSystem([:x => x_dynamics, :y => y_dynamics, :z => z_dynamics])
    test_dynamics(
        ds,
        IIP, D, M, DN, T,
        t0, vcat(SVector(x0), y0, z0), I,
        WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M))),
        @SMatrix [1.0 0.0 0.0 0.0 0.0 0.0 0.0
                  1.0 0.0 0.0 0.0 0.0 0.0 0.0
                  0.0 1.0 0.0 0.0 0.0 0.0 0.0
                  0.0 0.0 1.0 0.0 0.0 0.0 0.0
                  0.0 0.0 0.0 1.0 0.0 0.0 0.0
                  0.0 0.0 0.0 0.0 1.0 1.0 1.0
                  0.0 0.0 0.0 0.0 1.0 1.0 1.0
                  0.0 0.0 0.0 0.0 1.0 1.0 1.0
                  0.0 0.0 0.0 0.0 1.0 1.0 1.0]
    )
end
