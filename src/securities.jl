
abstract type AbstractSecurity end

# IDEA:
# me parece que llamarla SystemSecurity no esta mal pero me puede jugar en contra con el
# hecho que en los modelos tambien tengo que definir Security. Por otro lado, lo que podria
# hacer es que dentro de ModelSecurity se definan los SystemSecurities necesarios (pero por
# algo esto no lo hice antes...)
# El otro temita esta en el nombre, que es Security en singular pero puede representar
# objetos que son multidimensionales (aunque eso esta en System incluido). Por ahi habria
# que llamarlo SystemSecurities. Pero eso pronto lo veremos.
struct Security{Di,Df,Mi,Mf,U,X,DU,DX} <: AbstractSecurity
    u::U
    x::X
    du::DU
    dx::DX
    function Security{Di,Df,Mi,Mf}(u::U, x::X, du::DU, dx::DX) where {Di,Df,Mi,Mf,U,X,DU,DX}
        return new{Di,Df,Mi,Mf,U,X,DU,DX}(u, x, du, dx)
    end
end

function Security(
    ::AbstractDynamics{IIP,D,M}, d::Integer, m::Integer, u::AbstractVector, du::AbstractVector
) where {IIP,D,M}

    Di, Df = d, d + D - 1
    Mi, Mf = m, m + M - 1

    # la verdad es que en realidad no importan los views aca, ya que luego uso remake.
    # aca lo que me interesa es completar Di, Df, Mi y Mf
    x = view(u, Di:Df)
    dx = view(du, Di:Df)

    return Security{Di,Df,Mi,Mf}(u, x, du, dx)
end

function Security(
    ::AbstractDynamics{IIP,D,M}, d::Integer, m::Integer, u::AbstractVector, du::AbstractMatrix
) where {IIP,D,M}

    Di, Df = d, d + D - 1
    Mi, Mf = m, m + M - 1

    # la verdad es que en realidad no importan los views aca, ya que luego uso remake.
    # aca lo que me interesa es completar Di, Df, Mi y Mf
    x = view(u, Di:Df)
    dx = view(du, Di:Df, Mi:Mf)

    return Security{Di,Df,Mi,Mf}(u, x, du, dx)
end

function Security(
    ::AbstractDynamics{IIP,D,M}, d::Integer, m::Integer, u::AbstractVector, ::Nothing
) where {IIP,D,M}

    Di, Df = d, d + D - 1
    Mi, Mf = m, m + M - 1

    # la verdad es que en realidad no importan los views aca, ya que luego uso remake.
    # aca lo que me interesa es completar Di, Df, Mi y Mf
    x = view(u, Di:Df)

    return Security{Di,Df,Mi,Mf}(u, x, nothing, nothing)
end

function remake(::Security{Di,Df,Mi,Mf}, u::AbstractVector) where {Di,Df,Mi,Mf}
    return Security{Di,Df,Mi,Mf}(u, view(u, Di:Df), nothing, nothing)
end

function remake(::Security{Di,Df,Mi,Mf}, u::AbstractVector, du::AbstractVector) where {Di,Df,Mi,Mf}
    return Security{Di,Df,Mi,Mf}(u, view(u, Di:Df), du, view(du, Di:Df))
end

function remake(::Security{Di,Df,Mi,Mf}, u::AbstractVector, du::AbstractMatrix) where {Di,Df,Mi,Mf}
    return Security{Di,Df,Mi,Mf}(u, view(u, Di:Df), du, view(du, Di:Df, Mi:Mf))
end

(s::Security)() = s.x
(s::Security{D,D,M,M})() where {D,M} = s.x[]
(s::Security)(::Real) = s.x
(s::Security{D,D,M,M})(::Real) where {D,M} = s.x[]


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



abstract type ModelSecurity <: AbstractSecurity end

# include("model-security/equity.jl")

include("model-security/interest_rate.jl")
export FixedIncomeSecurities

# include("model-security/local_volatility.jl")
# include("model-security/stochastic_volatility.jl")