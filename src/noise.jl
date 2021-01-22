"""
    AbstractNoise{M}

Supertype for noise types with dimension `M`.
"""
abstract type AbstractNoise{M} end

dimension(::AbstractNoise{M}) where {M} = M

@doc raw"""
    ScalarNoise <: AbstractNoise{1}

A System of SDEs has `ScalarNoise` when there is only a unique noise process that affects
all the Stochastic Differential Equations. In that sense, the system is given by:
```math
dX⃗(t) = μ(X⃗(t), t) ⋅ dt + σ(X⃗(t), t) ⋅ dW(t),
```
with ``μ`` and ``σ`` functions that returns vectors... and ``dW`` a scalar function.
TODO: Obviamente escribir esto mejor.
"""
struct ScalarNoise <: AbstractNoise{1} end

@doc raw"""
    DiagonalNoise{M} <: AbstractNoise{M}

A system of SDEs:

```math
d\boldsymbol{x}(t) = d\boldsymbol{μ}(\boldsymbol{x}(t), t) \cdot dt + \dots,
```

has `DiagonalNoise` if ``\boldsymbol{σ}`` is a square diagonal matrix. (igualmente la
literatura asume que ademas,  la derivada σ[i,i] con respect a x_j es cero. sin embargo esto
solo cambia el orden de convergencia del scheme si es que no se cumple).

Aclarar que aqui simplificamos el asunto y usamos un sigma que devuelve un vector ya que es
lo mismo que usar la diagonal de una matriz.
TODO: Obviamente escribir esto mejor.
"""
struct DiagonalNoise{M} <: AbstractNoise{M} end

"""
    NonDiagonalNoise{M} <: AbstractNoise{M}

A System of SDEs has `NonDiagonalNoise`... tenemos que `D != M` o que σ no es una matriz
diagonal.
```math
dX⃗(t) = μ(X⃗(t), t) ⋅ dt + σ(X⃗(t), t) ⋅ dW⃗(t),
```
with ``μ`` and ``σ``
TODO: Obviamente escribir esto mejor.
"""
struct NonDiagonalNoise{M} <: AbstractNoise{M} end

DiagonalNoise(M::Integer) = DiagonalNoise{M}()
NonDiagonalNoise(M::Integer) = NonDiagonalNoise{M}()

function diffeqnoise(t0, ρ, IIP, D, M, DN)
    if DN
        if isequal(ρ, I)
            noise = nothing
        else
            if IIP
                W0 = zeros(M)
                Z0 = zeros(M)
                noise = CorrelatedWienerProcess!(ρ, t0, W0, Z0)
            else
                W0 = @SVector zeros(M)
                Z0 = @SVector zeros(M)
                noise = CorrelatedWienerProcess(ρ, t0, W0, Z0)
            end
        end
    else
        if isequal(D, M)
            if isequal(ρ, I)
                # NonDiagonalNoise with a square matrix and no correlations
                noise = nothing
            else
                if IIP
                    W0 = zeros(M)
                    Z0 = zeros(M)
                    noise = CorrelatedWienerProcess!(ρ, t0, W0, Z0)
                else
                    W0 = @SVector zeros(M)
                    Z0 = @SVector zeros(M)
                    noise = CorrelatedWienerProcess(ρ, t0, W0, Z0)
                end
            end
        else
            if isequal(ρ, I)
                # NonDiagonalNoise with non-square matrix and no correlations
                if IIP
                    W0 = zeros(M)
                    Z0 = zeros(M)
                    noise = WienerProcess!(t0, W0, Z0)
                else
                    W0 = @SVector zeros(M)
                    Z0 = @SVector zeros(M)
                    noise = WienerProcess(t0, W0, Z0)
                end
            else
                if IIP
                    W0 = zeros(M)
                    Z0 = zeros(M)
                    noise = CorrelatedWienerProcess!(ρ, t0, W0, Z0)
                else
                    W0 = @SVector zeros(M)
                    Z0 = @SVector zeros(M)
                    noise = CorrelatedWienerProcess(ρ, t0, W0, Z0)
                end
            end
        end
    end
    return noise
end