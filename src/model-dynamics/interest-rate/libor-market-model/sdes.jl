# TODO:
#  - HACER OOP VERSION Y VER SI EXPRESIONES VECTORIALES/MATRICIALES SON MAS RAPIDAS
#  - AGREGAR MEDIDAS CON NUMERAIRE P(t, Ti) para i arbitrario y fijo, i.e. Q^Ti measure
#  - HACER LAS CUENTAS MAS RAPIDAS CON LO DE RECURSIVIDAD
#  - CODEAR EL CASO NO DIAGONAL, ver que el paper copado usa eso en el ej. de la pg 32.
#  - Agregar la Tenor structure a las fechas de simulacion via tstops.

# IDEA: Las funciones de _searchsortedfirst o _searchsortedlast de DifferentialEquations.jl
# son mas rapidas que las de base y podrian ser de utilidad en esta funcion.
"""
    q(Tenors, t)

Returns the index in the tenor structure `Tenors` such that `Tenors[q(t)-1] ≤ t <
Tenors[q(t)]`, which corresponds to a left continuous interpretation.
"""
function q(Tenors, t)

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

@inline function drift!(du, u, p::LMMP{true,D,D,true,Terminal}, t) where {D}
    @unpack Tenors, τ, ρ, cache = p
    @unpack σ = cache
    p.σ(σ, t)

    @inbounds begin
        # 'i' ranges from 2 to D-1 because L₁ and LN have μ = 0
        for i in 2:D-1
            du[i] = zero(eltype(du))
            if t ≤ Tenors[i]
                for j in i+1:D
                    du[i] += (ρ[i,j] * τ[j] * σ[j] * u[j]) / (1 + τ[j] * u[j])
                end
                du[i] *= (-σ[i] * u[i])
            end
        end
    end

    return nothing
end

@inline function drift(u, p::LMMP{false,D,D,true,Terminal}, t) where {D}
    @unpack Tenors, τ, ρ = p
    σ = p.σ(t)

    du = @MVector zeros(eltype(u), D)

    @inbounds begin
        # 'i' ranges from 2 to D-1 because L₁ and LD have μ = 0
        for i in 2:D-1
            if t ≤ Tenors[i]
                for j in i+1:D
                    du[i] += (ρ[i,j] * τ[j] * σ[j] * u[j]) / (1 + τ[j] * u[j])
                end
                du[i] *= (-σ[i] * u[i])
            end
        end
    end

    return convert(SVector, m)
end

@inline function drift!(du, u, p::LMMP{true,D,D,true,Spot}, t) where {D}
    @unpack Tenors, τ, ρ, cache = p
    @unpack σ = cache
    p.σ(σ, t)

    # get left continuous index
    qt = q(Tenors, t)

    @inbounds begin
        # 'i' ranges from 2 to D-1 because L₁ and LD have μ = 0
        for i in 2:D-1
            du[i] = zero(eltype(du))
            if !isempty(qt:i)
                for j in qt:i
                    du[i] += (ρ[i,j] * τ[j] * σ[j] * u[j]) / (1 + τ[j] * u[j])
                end
                du[i] *= (σ[i] * u[i])
            end
        end
    end

    return nothing
end

@inline function drift(u, p::LMMP{false,D,D,true,Spot}, t) where {D}
    @unpack Tenors, τ, ρ = p
    σ = p.σ(t)

    du = @MVector zeros(eltype(u), D)

    # get left continuous index
    qt = q(Tenors, t)

    @inbounds begin
        # 'i' ranges from 2 to D-1 because L₁ and LD have μ = 0
        for i in 2:D-1
            if !isempty(qt:i)
                for j in qt:i
                    du[i] += (ρ[i,j] * τ[j] * σ[j] * u[j]) / (1 + τ[j] * u[j])
                end
                du[i] *= (σ[i] * u[i])
            end
        end
    end

    return convert(SVector, du)
end

@inline function diffusion!(du, u, p::LMMP{true,D,D,true}, t) where {D}
    @unpack Tenors, cache = p
    @unpack σ = cache
    p.σ(σ, t)

    @inbounds begin
        for i in 2:D
            du[i] = t > Tenors[i] ? zero(eltype(du)) : σ[i] * u[i]
        end
    end

    return nothing
end

@inline function diffusion(u, p::LMMP{false,D,D,true}, t) where {D}
    @unpack Tenors, τ, ρ = p
    σ = p.σ(t)

    du = @MVector zeros(eltype(u), D)

    @inbounds begin
        for i in 2:D
            du[i] = t > Tenors[i] ? zero(eltype(du)) : σ[i] * u[i]
        end
    end

    return convert(SVector, du)
end
