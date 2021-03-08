# TODO: add non-diagonal case from paper copado.

using UnPack

# setup from https://fbe.unimelb.edu.au/__data/assets/pdf_file/0009/2591829/200.pdf, section
# 8.2.

@testset "Terminal and Spot measures (IIP and OOP cases)" begin

for measure in (Terminal(), Spot())

    Δ = 1/2
    τ = @SVector [Δ for i in 1:12]
    Tenors = vcat(zero(eltype(τ)), cumsum(τ))

    # "abcd" time dependent volatility structure:
    function σt(i, t)
        return (0.05 + 0.09 * (Tenors[i] - t)) * exp(-0.44 * (Tenors[i] - t)) + 0.2
    end

    #! deberia ser una funcion de (p, t) y la version IIP de (u, p, t)?
    function σ(t)

        # this one computes unnecesary values
        # return @SVector [σt(i, t) for i in 1:4]

        # we could use an MVector, modify and convert to SVector

        # or this method:
        # notar que podria loopear solo en los indices 2:12-1 en Terminal measure
        return SVector(ntuple(Val{12}()) do i
            if t ≤ Tenors[i]
                return σt(i, t)
            else
                return zero(promote_type(1/t))
            end
        end)
    end

    ρ = @SMatrix [exp(-0.0669 * abs(i - j)) for i in 1:12, j in 1:12]
    L0 = @SVector [0.008 + 0.002 * i + 0.01 for i in 1:12]
    L = LiborMarketModelDynamics(L0, τ, σ, ρ, measure=measure, imethod=Schlogl(true))

    function f(u, p, t)
        @unpack _dynamics, _securities_ = p
        @unpack _L = _dynamics
        @unpack _L_ = _securities_

        L = remake(_L_, u, t)

        dL = drift(L(t), get_parameters(_L), t)

        return dL
    end

    function g(u, p, t)
        @unpack _dynamics, _securities_ = p
        @unpack _L = _dynamics
        @unpack _L_ = _securities_

        L = remake(_L_, u, t)

        dL = diffusion(L(t), get_parameters(_L), t)

        return dL
    end

    dynamics = [:L => L]
    ds_oop = DynamicalSystem(f, g, dynamics, nothing)
    sol_oop = solve(ds_oop, 6., seed=1)


    Δ = 1/2
    τ = [Δ for i in 1:12]
    Tenors = vcat(zero(eltype(τ)), cumsum(τ))

    function σt(i, t)
        return (0.05 + 0.09 * (Tenors[i] - t)) * exp(-0.44 * (Tenors[i] - t)) + 0.2
    end

    function σ!(u, t)

        # podria loopear solo en los indices 2:12-1 en Terminal measure
        for i in 1:12
            u[i] = zero(eltype(u))
            if t ≤ Tenors[i]
                u[i] = σt(i, t)
            end
        end

        return nothing
    end

    ρ = [exp(-0.0669 * abs(i - j)) for i in 1:12, j in 1:12]
    L0 = [0.008 + 0.002 * i + 0.01 for i in 1:12]
    L = LiborMarketModelDynamics(L0, τ, σ!, ρ, measure=measure, imethod=Schlogl(true))

    function f!(du, u, p, t)
        @unpack _dynamics, _securities_ = p
        @unpack _L = _dynamics
        @unpack _L_ = _securities_

        L = remake(_L_, du, u, t)

        drift!(L.dx, L(t), get_parameters(_L), t)

        return nothing
    end

    function g!(du, u, p, t)
        @unpack _dynamics, _securities_ = p
        @unpack _L = _dynamics
        @unpack _L_ = _securities_

        L = remake(_L_, du, u, t)

        diffusion!(L.dx, L(t), get_parameters(_L), t)

        return nothing
    end

    dynamics = [:L => L]
    ds_iip = DynamicalSystem(f!, g!, dynamics, nothing)
    sol_iip = solve(ds_iip, 6., seed=1)

    @test sol_oop.u ≈ sol_iip.u atol=1e-4
end
end
