
import DiffEqBase: remake

abstract type AbstractSecurity end

# FIXME: porque estoy guardando `u` y `du`?
struct Security{Di,Df,Mi,Mf,U,X,DU,DX} <: AbstractSecurity
    u::U
    x::X
    du::DU
    dx::DX
    function Security{Di,Df,Mi,Mf}(u::U, x::X, du::DU, dx::DX) where {Di,Df,Mi,Mf,U,X,DU,DX}
        return new{Di,Df,Mi,Mf,U,X,DU,DX}(u, x, du, dx)
    end
end

function Security(::AbstractDynamics{IIP,D,M}, d::Integer, m::Integer) where {IIP,D,M}
    Di, Df = d, d + D - 1
    Mi, Mf = m, m + M - 1
    return Security{Di,Df,Mi,Mf}(ntuple(_ -> nothing, 4)...)
end

function remake(::Security{Di,Df,Mi,Mf}, u::AbstractVector, t::Real) where {Di,Df,Mi,Mf}
    return Security{Di,Df,Mi,Mf}(
        u,
        s -> in_domain(s, t, view(u, Di:Df)),
        nothing,
        nothing
    )
end

function remake(::Security{Di,Df,Mi,Mf}, du::AbstractVector, u::AbstractVector, t::Real) where {Di,Df,Mi,Mf}
    return Security{Di,Df,Mi,Mf}(
        u,
        s -> in_domain(s, t, view(u, Di:Df)),
        du,
        view(du, Di:Df)
    )
end

function remake(::Security{Di,Df,Mi,Mf}, du::AbstractMatrix, u::AbstractVector, t::Real) where {Di,Df,Mi,Mf}
    return Security{Di,Df,Mi,Mf}(
        u,
        s -> in_domain(s, t, view(u, Di:Df)),
        du,
        view(du, Di:Df, Mi:Mf)
    )
end

(s::Security)(t::Real) = s.x(t)
(s::Security)(i::Integer, t::Real) = getindex(s.x(t), i)
(s::Security{D,D,M,M})(t::Real) where {D,M} = getindex(s.x(t), 1)

in_domain(s, t, result) = isequal(s, t) ? result : error("Outside time domain.")


#! IMPORTANTE:
# como ver la diagonal de una matrix?
# @btime @view (@view m[:])[diagind(m)]
# esto me va a servir cuando quiera escribir sobre el du de un dynamical system que no es DN
# pero que tenga un modelo de short rate que es DN
# Basicamente al crear la Security, lo que hago es comparar el DN del dynamical system con
# el DN del AbstractDynamics. Si son diferentes y estoy en IIP, tengo que hacer que el view
# sea sobre la diagonal de du.
# struct Security{DNAbastractDynamics,DNDynamicalSystem,...} end
# (s::Security{false,false,...})(du::AbstractMatrix) = view(du, Di:Df, Mi:Mf)
# (s::Security{true,false,...})(du::AbstractMatrix) = @view (@view du[:])[diagind(du)]
# donde en diagind(du) en realidad tengo que usar los indices que corresponden a esta AbstractDynamics
# Por ahi baste con hacer un remake del ShortRateModelDynamics (del AbstractDynamics en cuestion)



abstract type ModelSecurity <: AbstractSecurity end

# include("model-security/equity.jl")

include("model-security/interest_rate.jl")
export FixedIncomeSecurities

# include("model-security/local_volatility.jl")
# include("model-security/stochastic_volatility.jl")