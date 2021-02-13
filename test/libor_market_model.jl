using UnPack
using StochasticDiffEq

@testset "Terminal and Spot measures (IIP and OOP cases)" begin

for measure in (Terminal(), Spot())

    Δ = 0.25
    τ = @SVector [Δ for i in 1:4]
    Tenors = vcat(zero(eltype(τ)), cumsum(τ))

    # esto es de prueba viendo ej. pg 32 paper copado, aunque ahi es non-diagonal noise
    function σt(i, t)
        return 0.6 * exp(-0.1 * (Tenors[i] - t))
    end

    function σ(t)

        # this one computes unnecesary values
        # return @SVector [σt(i, t) for i in 1:4]

        # we could use an MVector, modify and convert to SVector

        # or this method:
        # notar que podria loopear solo en los indices 2:4-1 en Terminal measure
        return SVector(ntuple(Val{4}()) do i
            if t ≤ Tenors[i]
                return σt(i, t)
            else
                return zero(promote_type(1/t))
            end
        end)
    end

    ρ = @SMatrix [1.0 0.2 0.2 0.2
                  0.2 1.0 0.2 0.2
                  0.2 0.2 1.0 0.2
                  0.2 0.2 0.2 1.0]
    L0 = @SVector [0.0112, 0.0118, 0.0123, 0.0127]
    L = LiborMarketModelDynamics(L0, τ, σ, ρ, measure=measure, imethod=Schlogl(true))

    function f(u, p, t)
        @unpack L_dynamics, L_security = p

        L = remake(L_security, u)

        IR = FixedIncomeSecurities(L_dynamics, L)

        dL = UniversalDynamics.drift(L(t), UniversalDynamics.parameters(L_dynamics), t)

        return dL
    end

    function g(u, p, t)
        @unpack L_dynamics, L_security = p

        L = remake(L_security, u)

        IR = FixedIncomeSecurities(L_dynamics, L)

        dL = UniversalDynamics.diffusion(L(t), UniversalDynamics.parameters(L_dynamics), t)

        return dL
    end

    dynamics = OrderedDict(:L => L)
    ds_oop = DynamicalSystem(f, g, dynamics, nothing)
    sol_oop = solve(ds_oop, 1., seed=1)


    Δ = 0.25
    τ = [Δ for i in 1:4]
    Tenors = vcat(zero(eltype(τ)), cumsum(τ))

    # esto es de prueba viendo ej. pg 32 paper copado, aunque ahi es non-diagonal noise
    function σt(i, t)
        return 0.6 * exp(-0.1 * (Tenors[i] - t))
    end

    function σ!(u, t)

        # podria loopear solo en los indices 2:4-1 en Terminal measure
        for i in 1:4
            u[i] = zero(eltype(u))
            if t ≤ Tenors[i]
                u[i] = σt(i, t)
            end
        end

        return nothing
    end

    ρ = [1.0 0.2 0.2 0.2
        0.2 1.0 0.2 0.2
        0.2 0.2 1.0 0.2
        0.2 0.2 0.2 1.0]
    L0 = [0.0112, 0.0118, 0.0123, 0.0127]
    L = LiborMarketModelDynamics(L0, τ, σ!, ρ, measure=measure, imethod=Schlogl(true))

    function f!(du, u, p, t)
        @unpack L_dynamics, L_security = p

        L = remake(L_security, u, du)

        IR = FixedIncomeSecurities(L_dynamics, L)

        UniversalDynamics.drift!(L.dx, L(t), UniversalDynamics.parameters(L_dynamics), t)

        return nothing
    end

    function g!(du, u, p, t)
        @unpack L_dynamics, L_security = p

        L = remake(L_security, u, du)

        IR = FixedIncomeSecurities(L_dynamics, L)

        UniversalDynamics.diffusion!(L.dx, L(t), UniversalDynamics.parameters(L_dynamics), t)

        return nothing
    end

    dynamics = OrderedDict(:L => L)
    ds_iip = DynamicalSystem(f!, g!, dynamics, nothing)
    sol_iip = solve(ds_iip, 1., seed=1)

    @test sol_oop.u ≈ sol_iip.u atol = 1e-5
end
end
