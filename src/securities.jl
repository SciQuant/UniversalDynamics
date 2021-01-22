
abstract type AbstractSecurity end

#! la construccion de una security conviene que reciba como argumento un SystemDynamics?
#! Pregunto para saber si de ahi puedo sacar el IIP, D, M, DN.
#! Puede ser que el IIP, D, M, DN del SystemDynamics no coincidan con el que quiero en
#! la security porque la security esta embebida dentro de un DynamicalSystem. Es decir,
#! IIP, D y M vienen del SystemDynamics, eso seguro. Pero DN viene del DynamicalSystem

# IDEA:
# me parece que llamarla SystemSecurity no esta mal pero me puede jugar en contra con el
# hecho que en los modelos tambien tengo que definir Security. Por otro lado, lo que podria
# hacer es que dentro de ModelSecurity se definan los SystemSecurities necesarios (pero por
# algo esto no lo hice antes...)
# El otro temita esta en el nombre, que es Security en singular pero puede representar
# objetos que son multidimensionales (aunque eso esta en System incluido). Por ahi habria
# que llamarlo SystemSecurities. Pero eso pronto lo veremos.
#! struct SystemSecurity{IIP,D,M,DN,S,T} <: AbstractSecurity
struct Security{IIP,D,M,DN,S,T} <: AbstractSecurity
    x::S
    dx::T
end

# Cases:
#   Security inside system drift with:
#     - Out Of Place
#     - OneDimensional case
#     - Noise type is irrelevant in drift, i.e. for any `M` and `DN`
function Security{false,1,M,DN}(u::SVector, idx::Integer, t::Real) where {M,DN}
    # xₜ = view(u, idx) # si uso view y desreferencio el functor, este metodo termina igual al de D > 1
    xₜ = u[idx]
    x = s -> isequal(s, t) ? xₜ : throw(DomainError(s, "time must be $(string(t)) instead of $(string(s))"))
    return Security{false,1,M,DN,typeof(x),Nothing}(x, nothing)
end

# same case as above but I believe que esta la necesito para ModelMacro
function Security{false,1,M,DN}(u::SVector, idxs::UnitRange, t::Real) where {M,DN}
    idxs.start == idxs.stop || throw(DomainError(idxs, "inconsistent indexes for a one dimensional security."))
    return Security{false,1,M,DN}(u, idxs.start, t)
end

# Cases:
#   Security inside system drift with:
#     - Out Of Place
#     - MultiDimensional dimensional case
#     - Noise type is irrelevant in drift, i.e. for any `M` and `DN`
function Security{false,D,M,DN}(u::SVector, idxs::UnitRange, t::Real) where {D,M,DN}
    xₜ = view(u, idxs)
    x = s -> isequal(s, t) ? xₜ : throw(DomainError(s, "time must be $(string(t)) instead of $(string(s))"))
    return Security{false,D,M,DN,typeof(x),Nothing}(x, nothing)
end

(s::Security{false})(t::Real) = s.x(t)

#! TO BE COMPLETED: faltan banda de casos que ire haciendo.


abstract type ModelSecurity <: AbstractSecurity end

# include("model-security/equity.jl")

include("model-security/interest_rate.jl")
export FixedIncomeSecurities

# include("model-security/local_volatility.jl")
# include("model-security/stochastic_volatility.jl")