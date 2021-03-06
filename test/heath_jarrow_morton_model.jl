using UnPack


@testset "Terminal and Spot measures (IIP and OOP cases)" begin
for measure in (Terminal(), Spot())
    Δ = 0.25
    τ = @SVector [Δ for i in 1:4]
    Tenors = vcat(zero(eltype(τ)), cumsum(τ))

    α₁ = 0.10
    α₂ = 0.13

    σ₁₁ = 0.11
    σ₂₁ = 0.12
    σ₂₂ = 0.14

    # Volatility as functions
    σ₁(t, T) = σ₁₁ * exp(- α₁ * (T - t))
    σ₂(t, T) = σ₂₁ * exp(- α₁ * (T - t)) + σ₂₂ * exp(- α₂ * (T - t))

    σt(t, T) = @SVector [σ₁(t, T), σ₂(t, T), σ₁(t, T), σ₁(t, T)]

    function σ(t, T)
        return SVector(ntuple(Val{4}()) do i
            if t ≤ T
                # porque no mandas `i` como argumento de la funcion y calculas en base a
                # eso en lugar de calcular todos y tomar la componente `i`? lo podes hacer
                # definiendo los vectores/matrices σ y α correspondientes... y si necesitas
                # una componente `j` para hacerlo mejor aun, hacelo...
                return σt(t, T)[i]
            else
                return zero(promote_type(1/t))
            end
        end)
    end


    F0 = @SVector [0.0112, 0.0118, 0.0123, 0.0127]
    F = HeathJarrowMortonModelDynamics(F0, τ, σ, measure=measure)

    function f(u, p, t)
        @unpack _dynamics, _securities_ = p
        @unpack _F = _dynamics
        @unpack _F_ = _securities_

        F = remake(_F_, u, t)

        dF = drift(F(t), get_parameters(_F), t)

        return dF
    end

    function g(u, p, t)
        @unpack _dynamics, _securities_ = p
        @unpack _F = _dynamics
        @unpack _F_ = _securities_

        F = remake(_F_, u, t)

        dF = diffusion(F(t), get_parameters(_F), t)

        return dF
    end

    dynamics = [:F => F]
    ds_oop = DynamicalSystem(f, g, dynamics, nothing)
    sol_oop = solve(ds_oop, 1., seed=1)

    function σ!(u, t, T)
        for i in 1:4
            u[i] = zero(eltype(u))
            if t ≤ T
                #! FIXME: idem al comentario de arriba...
                u[i] = σt(t, T)[i]
            end
        end

        return nothing
    end

    F0 = [0.0112, 0.0118, 0.0123, 0.0127]
    F1 = HeathJarrowMortonModelDynamics(F0, τ, σ!, measure=measure)

    function f!(du, u, p, t)
        @unpack _dynamics, _securities_ = p
        @unpack _F = _dynamics
        @unpack _F_ = _securities_

        F = remake(_F_, du, u, t)

        drift!(F.dx, F(t), get_parameters(_F), t)

        return nothing
    end

    function g!(du, u, p, t)
        @unpack _dynamics, _securities_ = p
        @unpack _F = _dynamics
        @unpack _F_ = _securities_

        F = remake(_F_, du, u, t)

        diffusion!(F.dx, F(t), get_parameters(_F), t)

        return nothing
    end

    dynamics = [:F => F1]
    ds_iip = DynamicalSystem(f!, g!, dynamics, nothing)
    sol_iip = solve(ds_iip, 1., seed=1)

    @test sol_oop.u ≈ sol_iip.u atol=1e-4

end
end