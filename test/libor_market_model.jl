
@testset "Terminal measure" begin

    N = 4
    Δ = 0.25
    τ = @SVector [Δ for i in 1:N]
    Tenors = vcat(zero(eltype(τ)), cumsum(τ))

    # esto es de prueba viendo ej. pg 32 paper copado, aunque ahi es non-diagonal noise
    function σt(i, t)
        return 0.6 * exp(0.1 * (Tenors[i] - t))
    end

    function σ(t)

        # this one computes unnecesary values
        # return @SVector [σt(i, t) for i in 1:4]

        # use MVector, modify and convert to SVector

        # or this method:
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
    L = LiborMarketModelDynamics(L0, τ, σ, ρ, measure=Terminal(), imethod=Schlogl(true))

    function f(u, p, t)
        @unpack L_dynamics, L_security = p

        L = UniversalDynamics.remake(L_security, u)

        IR = FixedIncomeSecurities(L_dynamics, L)

        dL = UniversalDynamics.drift(L(t), UniversalDynamics.parameters(L_dynamics), t)

        return dL
    end

    function g(u, p, t)
        @unpack L_dynamics, L_security = p

        L = UniversalDynamics.remake(L_security, u)

        IR = FixedIncomeSecurities(L_dynamics, L)

        dL = UniversalDynamics.drift(L(t), UniversalDynamics.parameters(L_dynamics), t)

        return dL
    end

    dynamics = OrderedDict(:L => L)
    ds_oop = DynamicalSystem(f, g, dynamics, nothing)

    # sol_oop = solve(ds_oop, 1., seed=1) # incluso cambiando el ruido para que tenga w0 y z0 con Array aun no lo pude solucionar
end