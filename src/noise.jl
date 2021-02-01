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
du⃗(t) = f(t, u⃗(t)) ⋅ dt + g(t, u⃗(t)) ⋅ dW(t), \quad u⃗(t₀) = u⃗₀,
```
in a time span ``\mathbb{I} = \left[ t₀, T \right]`` with drift coefficient ``f \colon
\mathbb{I} × \mathbb{R}ᴰ → \mathbb{R}ᴰ``, diffusion coefficient ``g \colon \mathbb{I} ×
\mathbb{R}ᴰ → \mathbb{R}ᴰ``, ``1``-dimensional driving Wiener process ``dW(t)`` and initial
condition ``u⃗₀``.
"""
struct ScalarNoise <: AbstractNoise{1} end

@doc raw"""
    DiagonalNoise{M} <: AbstractNoise{M}

A ``D``-dimensional system has `DiagonalNoise` when there are ``M = D`` noise processes that
affect each Stochastic Differential Equation individually. In that sense, the system is
given by:
```math
du⃗(t) = f(t, u⃗(t)) ⋅ dt + g(t, u⃗(t)) ⋅ d\vec{W}(t), \quad u⃗(t₀) = u⃗₀,
```
in a time span ``\mathbb{I} = \left[ t₀, T \right]`` with drift coefficient ``f \colon
\mathbb{I} × \mathbb{R}ᴰ → \mathbb{R}ᴰ``, diffusion coefficient ``g \colon \mathbb{I} ×
\mathbb{R}^D → \mathrm{diag} \colon \mathbb{R}^{D × D}``, ``D``-dimensional driving Wiener
correlated or uncorrelated process ``d\vec{W}(t)`` and initial condition ``u⃗₀``.

In these kind of systems, the diagonal of ``g(t, u⃗(t))`` is represented by a
``D``-dimensional vector. Then, the product ``g(t, u⃗(t)) ⋅ d\vec{W}(t)`` is replaced by the
broadcasted product `.*`.
"""
struct DiagonalNoise{M} <: AbstractNoise{M} end

DiagonalNoise{1}() = ScalarNoise()

DiagonalNoise(M::Integer) = DiagonalNoise{M}()

@doc raw"""
    NonDiagonalNoise{M} <: AbstractNoise{M}

A ``D``-dimensional system has `NonDiagonalNoise` when there are M noise processes that
affect the Stochastic Differential Equations. In that sense, the system is given by:
```math
du⃗(t) = f(t, u⃗(t)) ⋅ dt + g(t, u⃗(t)) ⋅ d\vec{W}(t), \quad u⃗(t₀) = u⃗₀,
```
in a time span ``\mathbb{I} = \left[ t₀, T \right]`` with drift coefficient ``f \colon
\mathbb{I} × \mathbb{R}ᴰ → \mathbb{R}ᴰ``, diffusion coefficient ``g \colon \mathbb{I} ×
\mathbb{R}ᴰ → \mathbb{R}^{D × M}``, ``M``-dimensional driving Wiener correlated or
uncorrelated process ``d\vec{W}(t)`` and initial condition ``u⃗₀``.
"""
struct NonDiagonalNoise{M} <: AbstractNoise{M} end

NonDiagonalNoise{1}() = ScalarNoise()

NonDiagonalNoise(M::Integer) = NonDiagonalNoise{M}()


function diffeqnoise(t0, ρ, IIP, D, M, DN, ep=true)
    if DN
        if isequal(ρ, I)
            noise = nothing
        else
            if IIP
                W0 = zeros(M)
                Z0 = ep ? zeros(M) : nothing
                noise = CorrelatedWienerProcess!(ρ, t0, W0, Z0)
            else
                W0 = @SVector(zeros(M))
                Z0 = ep ? @SVector(zeros(M)) : nothing
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
                    Z0 = ep ? zeros(M) : nothing
                    noise = CorrelatedWienerProcess!(ρ, t0, W0, Z0)
                else
                    W0 = @SVector(zeros(M))
                    Z0 = ep ? @SVector(zeros(M)) : nothing
                    noise = CorrelatedWienerProcess(ρ, t0, W0, Z0)
                end
            end
        else
            if isequal(ρ, I)
                # NonDiagonalNoise with non-square matrix and no correlations
                if IIP
                    W0 = zeros(M)
                    Z0 = ep ? zeros(M) : nothing
                    noise = WienerProcess!(t0, W0, Z0)
                else
                    W0 = @SVector(zeros(M))
                    Z0 = ep ? @SVector(zeros(M)) : nothing
                    noise = WienerProcess(t0, W0, Z0)
                end
            else
                if IIP
                    W0 = zeros(M)
                    Z0 = ep ? zeros(M) : nothing
                    noise = CorrelatedWienerProcess!(ρ, t0, W0, Z0)
                else
                    W0 = @SVector(zeros(M))
                    Z0 = ep ? @SVector(zeros(M)) : nothing
                    noise = CorrelatedWienerProcess(ρ, t0, W0, Z0)
                end
            end
        end
    end
    return noise
end

function diffeq_noise_rate_prototype(IIP, D, M, DN, dynamics)
    if DN
        # DiagonalNoise
        noise_rate_prototype = nothing

    elseif isone(M)
        # ScalarNoise:
        # El noise rate prototype sera una matrix de una columna y D filas (i.e. un vector
        # de dimension D). Sin embargo, en este caso el noise_rate_prototype sale de la
        # condicion inicial u0, por lo que no hace falta enviar noise_rate_prototype.
        noise_rate_prototype = nothing

    else
        prototypes = []
        for dynamic in dynamics
            if diagonalnoise(dynamic)
                # For DiagonalNoise, gprototype is a vector. But, when mixed with other
                # dynamics, it should be casted to a diagonal matrix.
                push!(prototypes, Diagonal(gprototype(dynamic)))
            else
                push!(prototypes, gprototype(dynamic))
            end
        end
        noise_rate_prototype = cat(prototypes...; dims=(1, 2))

        if IIP
            # We cast to Matrix because cat of a Diagonal matrix returns a sparse object.
            # However, in some cases we might want to return sparse matrices here. We must
            # allow that.
            noise_rate_prototype = convert(Matrix,noise_rate_prototype)
        else
            noise_rate_prototype = convert(SMatrix{D,M}, noise_rate_prototype)
        end
    end

    return noise_rate_prototype
end