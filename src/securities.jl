
abstract type AbstractSecurity end

# IDEA:
# me parece que llamarla SystemSecurity no esta mal pero me puede jugar en contra con el
# hecho que en los modelos tambien tengo que definir Security. Por otro lado, lo que podria
# hacer es que dentro de ModelSecurity se definan los SystemSecurities necesarios (pero por
# algo esto no lo hice antes...)
# El otro temita esta en el nombre, que es Security en singular pero puede representar
# objetos que son multidimensionales (aunque eso esta en System incluido). Por ahi habria
# que llamarlo SystemSecurities. Pero eso pronto lo veremos.
struct Security{Di,Df,Mi,Mf} <: AbstractSecurity end

function Security(sd::AbstractDynamics, d::Integer, m::Integer)
    D = dimension(sd)
    M = noise_dimension(sd)

    Di, Df = d, d + D - 1
    Mi, Mf = m, m + M - 1

    return Security{Di,Df,Mi,Mf}()
end

(s::Security{Di,Df})(u::AbstractVector) where {Di,Df} = view(u, Di:Df)
(s::Security{Di,Df,Mi,Mf})(du::AbstractMatrix) where {Di,Df,Mi,Mf} = view(du, Di:Df, Mi:Mf)

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