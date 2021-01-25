# TODO:
#  - HACER OOP VERSION Y VER SI EXPRESIONES VECTORIALES/MATRICIALES SON MAS RAPIDAS
#  - AGREGAR MEDIDAS CON NUMERAIRE P(t, Ti) para i arbitrario y fijo, i.e. Q^Ti measure
#  - HACER LAS CUENTAS MAS RAPIDAS CON LO DE RECURSIVIDAD
#  - CODEAR EL CASO NO DIAGONAL, ver que el paper copado usa eso en el ej. de la pg 32.
#  - Agregar la Tenor structure a las fechas de simulacion via tstops.

# IDEA: Las funciones de _searchsortedfirst o _searchsortedlast de DifferentialEquations.jl
# son mas rapidas que las de base y podrian ser de utilidad en esta funcion.
"""
    q(Tenors::Vector{<:Real}, t::Real)

Returns the index in the tenor structure `Tenors` such that `Tenors[q(t)-1] ≤ t <
Tenors[q(t)]`, which corresponds to a left continuous interpretation.
"""
function q(Tenors::AbstractVector{<:Real}, t::Real)

    if t > Tenors[end] || t < Tenors[begin]
        throw(DomainError("`t` must lie in `Tenors`."))
    end

    # puede caer justo sobre un elemento o dentro de un intevalo
    r = searchsorted(Tenors, t)

    if !isempty(r)
        # nos fijamos si es el ultimo porque a ese lo tratamos distinto
        return isequal(t, Tenors[end]) ? getindex(r) : getindex(r) + 1
    else
        # si cae sobre un intervalo, devolvemos el indice por derecha
        return first(r)
    end
end

#! IMPORTARNTE
# TODO: los dispatchs podrian ser mas prolijos, por ejemplo:
# function drift!(du, u, p, t)
#     IIP = isinplace(p)
#     D = dimension(p)
#     M = noise_dimension(p)
#     DN = diagonalnoise(p)
#     M = measure(p)
#     # S = # sigma type
#     return drift!(Val{IIP}, ...)
# end
#! Los dispatchs como estan ahora estan deprecarted porque cambie el orden de las cosas, asi
#! que no van a funcionar (o si funcionan estan mal). Hay que hacerlo con lo de arriba si o
#! si.

@inline function drift!(
    du::AbstractVector{<:Real},
    u::AbstractVector{<:Real},
    p::LMMP{true,true,Terminal,<:AbstractVector{<:Real}},
    t::Real
)
    @unpack Tenors, τ, σ, ρ, N = p

    z = zero(eltype(du))

    @inbounds begin
        # 'i' ranges from 2 to N-1 because L₁ and LN have μ = 0
        for i = 2:N-1
            res = z

            if t > Tenors[i]
                du[i] = res
            else
                # TODO: investigate if we can avoid allocations in sum()
                for j in i+1:N
                    res += (ρ[i, j] * τ[j] * σ[j] * u[j]) / (1 + τ[j] * u[j])
                end
                du[i] = -σ[i] * u[i] * res
            end
        end
    end

    return nothing
end

@inline function drift!(
    du::AbstractVector{<:Real},
    u::AbstractVector{<:Real},
    p::LMMP{true,true,Terminal,F},
    t::Real
) where {F<:Function}
    @unpack Tenors, τ, σ, ρ, N = p

    z = zero(eltype(du))

    @inbounds begin
        # 'i' ranges from 2 to N-1 because L₁ and LN have μ = 0
        for i in 2:N-1
            res = z

            if t > Tenors[i]
                du[i] = res
            else
                # TODO: investigate if we can avoid allocations in sum()
                for j in i+1:N
                    res += (ρ[i, j] * τ[j] * σ(j, t) * u[j]) / (1 + τ[j] * u[j])
                end
                du[i] = -σ(i, t) * u[i] * res
            end
        end
    end

    return nothing
end

@inline function drift!(
    du::AbstractVector{<:Real},
    u::AbstractVector{<:Real},
    p::LMMP{true,true,Spot,<:AbstractVector{<:Real}},
    t::Real
)
    @unpack Tenors, τ, σ, ρ, N = p

    z = zero(eltype(du))

    # get left continuous index
    qt = q(Tenors, t)

    @inbounds begin
        # 'i' ranges from 2 because L₁ have μ = 0
        for i in 2:N
            res = z
            if isempty(qt:i)
                du[i] = res
            else
                # TODO: investigate if we can avoid allocations in sum()
                for j in qt:i
                    res += (ρ[i, j] * τ[j] * σ[j] * u[j]) / (1 + τ[j] * u[j])
                end
                du[i] = σ[i] * u[i] * res
            end
        end
    end

    return nothing
end

@inline function drift!(
    du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, p::LMMP{true,true,Spot,F}, t::Real
) where {F<:Function}
    @unpack Tenors, τ, σ, ρ, N = p

    z = zero(eltype(du))

    qt = q(Tenors, t)

    @inbounds begin
        # 'i' ranges from 2 because L₁ have μ = 0
        for i in 2:N
            res = z
            if isempty(qt:i)
                du[i] = res
            else
                # TODO: investigate if we can avoid allocations in sum()
                for j in qt:i
                    res += (ρ[i, j] * τ[j] * σ(j, t) * u[j]) / (1 + τ[j] * u[j])
                end
                du[i] = σ(i, t) * u[i] * res
            end
        end
    end

    return nothing
end

@inline function diffusion!(
    du::AbstractVector{<:Real},
    u::AbstractVector{<:Real},
    p::LMMP{true,true,M,<:AbstractVector{<:Real}},
    t::Real
) where {M}
    @unpack Tenors, σ, N = p

    z = zero(eltype(du))

    @inbounds begin
        for i in 2:N
            du[i] = t > Tenors[i] ? z : σ[i] * u[i]
        end
    end

    return nothing
end

@inline function diffusion!(
    du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, p::LMMP{true,true,M,F}, t::Real
) where {M,F<:Function}
    @unpack Tenors, σ, N = p

    z = zero(eltype(du))

    @inbounds begin
        for i in 2:N
            du[i] = t > Tenors[i] ? z : σ(i, t) * u[i]
        end
    end

    return nothing
end