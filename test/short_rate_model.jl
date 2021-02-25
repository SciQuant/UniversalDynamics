using UnPack
using StochasticDiffEq

include("DaiSingletonParameters_A3_1.jl")

@testset "DaiSingleton 𝔸₁(3)ₘₐₓ" begin

    (υ₀, θ₀, r₀, μ, ν, κ_rυ, κ, ῡ, θ̄, η, σ_θυ, σ_θr, σ_rυ, σ_rθ, ζ, α_r, β_θ) = DaiSingletonParameters()

    x0 = @SVector [υ₀, θ₀, r₀]

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
        1   0 0
        β_θ 0 0
        1   0 0
    ]

    x = MultiFactorAffineModelDynamics(x0, ϰ, θ, Σ, α, β, ξ₀, ξ₁; noise=NonDiagonalNoise(3))
    B = SystemDynamics(one(eltype(x)))

    function f(u, p, t)
        @unpack _dynamics, _securities_ = p
        @unpack _x = _dynamics
        @unpack _x_, _B_ = _securities_

        x = remake(_x_, u)
        B = remake(_B_, u)

        IR = FixedIncomeSecurities(_x, x, B)

        dx = drift(x(t), get_parameters(_x), t)
        dB = IR.r(t) * B(t)

        return vcat(dx, dB)
    end

    function g(u, p, t)
        @unpack _dynamics, _securities_ = p
        @unpack _x = _dynamics
        @unpack _x_, _B_ = _securities_

        x = remake(_x_, u)
        B = remake(_B_, u)

        dx = diffusion(x(t), get_parameters(_x), t)
        dB = zero(eltype(u)) # @SMatrix zeros(eltype(u), 1, 1)

        return @SMatrix [dx[1,1] dx[1,2] dx[1,3]  0
                         dx[2,1] dx[2,2] dx[2,3]  0
                         dx[3,1] dx[3,2] dx[3,3]  0
                               0       0       0 dB]
    end

    dynamics = [:x => x, :B => B]
    ds_oop = DynamicalSystem(f, g, dynamics, nothing)
    sol_oop = solve(ds_oop, 1., alg=EM(), dt=0.01, seed=1)

    # para hacer funcionar esto tengo que sacar a f fuera del bloque de test y ahi veo que
    # no hay allocations
    # u = @SVector ones(4)
    # p = ds_oop.params
    # t = 0.05
    # @btime f($u, $p, $t)


    x0 = [υ₀, θ₀, r₀]

    ξ₀!(t) = zero(t) # ξ₀ = zero

    function ξ₁!(u, t)
        u[1] = 0
        u[2] = 0
        u[3] = 1
        return nothing
    end

    function ϰ!(u, t)
        u[1,1] = μ
        u[2,2] = ν
        u[3,1] = κ_rυ
        u[3,2] = -κ
        u[3,3] = κ
        return nothing
    end

    function θ!(u, t)
        u[1] = ῡ
        u[2] = θ̄
        u[3] = θ̄
        return nothing
    end

    function Σ!(u, t)
        u[1,1] = η
        u[2,1] = η * σ_θυ
        u[2,2] = 1
        u[2,3] = σ_θr
        u[3,1] = η * σ_rυ
        u[3,2] = σ_rθ
        u[3,3] = 1
        return nothing
    end

    function α!(u, t)
        u[1] = 0
        u[2] = ζ^2
        u[3] = α_r
        return nothing
    end

    function β!(u, t)
        u[1,1] = 1
        u[2,1] = β_θ
        u[3,1] = 1
        return nothing
    end

    # declare short rate model dynamics
    x = MultiFactorAffineModelDynamics(x0, ϰ!, θ!, Σ!, α!, β!, ξ₀!, ξ₁!; noise=NonDiagonalNoise(3))

    # declare money market account dynamics
    B = SystemDynamics(ones(eltype(x), 1))

    # in place drift coefficient
    function f!(du, u, p, t)
        @unpack _dynamics, _securities_ = p
        @unpack _x = _dynamics
        @unpack _x_, _B_ = _securities_

        x = remake(_x_, u, du)
        B = remake(_B_, u, du)

        IR = FixedIncomeSecurities(_x, x, B)

        drift!(x.dx, x(t), get_parameters(_x), t)
        B.dx[] = IR.r(t) * B(t)

        return nothing
    end

    # in place diffusion coefficient
    function g!(du, u, p, t)
        @unpack _dynamics, _securities_ = p
        @unpack _x = _dynamics
        @unpack _x_, _B_ = _securities_

        x = remake(_x_, u, du)
        B = remake(_B_, u, du)

        diffusion!(x.dx, x(t), get_parameters(_x), t)
        B.dx[] = zero(eltype(u))

        return nothing
    end

    dynamics = [:x => x, :B => B]
    ds_iip = DynamicalSystem(f!, g!, dynamics, nothing)
    sol_iip = solve(ds_iip, 1., alg=EM(), dt=0.01, seed=1)

    # idem con la de arriba
    # du = zeros(4)
    # u = ones(4)
    # p = ds_iip.params
    # t = 0.05
    # @btime $f!($du, $u, $p, $t)

    @test sol_oop.u ≈ sol_iip.u
end



