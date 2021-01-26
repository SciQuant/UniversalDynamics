"""
    AbstractNoise{M}

Supertype for noise types with dimension `M`.
"""
abstract type AbstractNoise{M} end

dimension(::AbstractNoise{M}) where {M} = M

@doc raw"""
    ScalarNoise <: AbstractNoise{1}

A ``D``-dimensional system has `ScalarNoise` when there is only a unique noise process that
affects all the Stochastic Differential Equations. In that sense, the system is given by:
```math
du⃗(t) = f(t, u⃗(t)) ⋅ dt + g(t, u⃗(t)) ⋅ dW(t), \quad u⃗(t_0) = u⃗₀,
```
with drift coefficient ``f \colon \left[ t₀, T \right] × \mathbb{R}ᴰ → \mathbb{R}ᴰ``,
diffusion coefficient ``g \colon \left[ t₀, T \right] × \mathbb{R}ᴰ → \mathbb{R}ᴰ``,
``1``-dimensional driving Wiener process ``dW(t)`` and initial condition ``u⃗₀``.
"""
struct ScalarNoise <: AbstractNoise{1} end

@doc raw"""
    DiagonalNoise{M} <: AbstractNoise{M}

A ``D``-dimensional system has `DiagonalNoise` when there are ``M = D`` noise processes that
affect each Stochastic Differential Equation individually. In that sense, the system is
given by:
```math
du⃗(t) = f(t, u⃗(t)) ⋅ dt + g(t, u⃗(t)) ⋅ d\vec{W}(t), \quad u⃗(t_0) = u⃗₀,
```
with drift coefficient ``f \colon \left[ t₀, T \right] × \mathbb{R}ᴰ → \mathbb{R}ᴰ``,
diffusion coefficient ``g \colon \left[ t₀, T \right] × \mathbb{R}^D → \mathrm{diag} \colon
\mathbb{R}^{D × D}``, ``D``-dimensional driving Wiener correlated or uncorrelated process
``d\vec{W}(t)`` and initial condition ``u⃗₀``.
"""
struct DiagonalNoise{M} <: AbstractNoise{M} end

DiagonalNoise{1}() = ScalarNoise()

DiagonalNoise(M::Integer) = DiagonalNoise{M}()

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

NonDiagonalNoise{1}() = ScalarNoise()

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