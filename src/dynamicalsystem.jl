"""
    DynamicalSystem{IIP,D,M,DN,T} <: AbstractDynamics{IIP,D,M,DN,T}

The central structure of **UniversalDynamics**. Its main purpose is to construct an
[`AbstractDynamics`](@ref) given a collection of [`AbstractDynamics`](@ref).

## Type parameters:
See [`AbstractDynamics`](@ref) for detailed information.

## Fields:
- `f`: drift coefficient represented as either an in place or out of place function,
- `g`: diffusion coefficient represented as either an in place or out of place function,
- `attributes`: see [`DynamicsAttributes`](@ref) for detailed information,
- `params`: container with parameters values (preferably a `NamedTuple`),
- `dynamics`: collection of [`AbstractDynamics`](@ref) (preferably in a `OrderedDict`),
- `securities`: useful handlers for both dynamics simulation and derivatives pricing.

## Declaration:

Given a collection of dynamics, group them together in a container such as an `OrderedDict`
and define a `DynamicalSystem` using:

```julia
DynamicalSystem(dynamics)
```

The resulting object will have information that can be inspected and is useful for coding
the coefficients (drift and difussion) functions.

The complete declaration of a `DynamicalSystem` requires:

```julia
DynamicalSystem(f, g, dynamics, params)
```

with `f` and `g` either the in place or the out of place functions for the coefficients:
- **Out of place**: coefficients must be in the form `f(u, p, t) -> SVector` and `g(u, p, t)
  -> SVector` for [`DiagonalNoise`](@ref) cases or `g(u, p, t) -> SMatrix` for
  [`NonDiagonalNoise`](@ref) cases. These functions return the drift or diffusion
  coefficients as `SArrays` given a current state `u`, current time `t` and a set of
  parameters.
- **In place**: coefficients must be in the form `f(du, u, p, t) -> nothing` and `g(du, u,
  p, t) -> nothing`. These functions modify in place `du::Array` and set it equal to either
  the drift or the diffusion coefficients given a current state `u`, current time `t` and a
  set of parameters.
"""
struct DynamicalSystem{IIP,D,M,DN,T,F,G,A,P,DS,S} <: AbstractDynamics{IIP,D,M,DN,T}
    f::F
    g::G
    attributes::A
    params::P
    dynamics::DS
    securities::S

    function DynamicalSystem{IIP,D,M,DN,T}(
        f::F, g::G, attrs::A, params::P, dynamics::DS, securities::S
    ) where {IIP,D,M,DN,T,F,G,A,P,DS,S}
        return new{IIP,D,M,DN,T,F,G,A,P,DS,S}(f, g, attrs, params, dynamics, securities)
    end
end

function DynamicalSystem(f, g, dynamics_container, params)

    if isempty(dynamics_container)
        throw(ArgumentError("provide at least one `AbstractDynamics`."))
    end

    dynamics = values(dynamics_container)

    IIP = all(isinplace.(dynamics))
    OOP = all((!isinplace).(dynamics))

    if isequal(IIP, OOP)
        error("all dynamics *must* be either In-Place or Out-Of-Place.")
    end

    D = sum(dimension.(dynamics))
    M = sum(noise_dimension.(dynamics))

    DN = !any((!diagonalnoise).(dynamics))

    t0s = [initialtime(d) for d in dynamics]
    if !all(t -> t == first(t0s), t0s)
        error("all dynamics *must* have the same initial time.")
    end
    t0 = initialtime(first(dynamics))

    T = typeof(t0)

    x0 = IIP ? vcat(state.(dynamics)...) : vcat(SVector.(state.(dynamics))...)

    ρ = cat(cor.(dynamics)..., dims = (1, 2))
    ρ = IIP ? Array{T,2}(ρ) : SMatrix{size(ρ)...,T}(ρ)

    noise_rate_prototype = diffeq_noise_rate_prototype(IIP, D, M, DN, dynamics)

    #! por ahi deberiamos definir el noise recien en la llamada a solve
    noise = diffeqnoise(t0, ρ, IIP, D, M, DN)

    attrs = DynamicsAttributes(t0, x0, ρ, noise, noise_rate_prototype)

    securities = Dict()
    d = m = 1
    for (name, abstract_dynamics) in dynamics_container
        x = Security(abstract_dynamics, d, m)
        push!(securities, Symbol(name, :_security) => x)
        d += dimension(abstract_dynamics)
        m += noise_dimension(abstract_dynamics)
    end

    securities = (; securities...)

    dynamics = Dict(Symbol(key, :_dynamics) => value for (key, value) in dynamics_container)
    dynamics = (; dynamics...)

    return DynamicalSystem{IIP,D,M,DN,T}(f, g, attrs, params, dynamics, securities)
end

DynamicalSystem(dynamics) = DynamicalSystem(nothing, nothing, dynamics, nothing)

for method in (:initialtime, :state, :cor, :noise, :noise_rate_prototype)
    @eval begin
        $method(ds::DynamicalSystem) = $method(ds.attributes)
    end
end

parameters(ds::DynamicalSystem) = parameters(ds, ds.params)
parameters(ds::DynamicalSystem, params) = merge(params, ds.dynamics, ds.securities)
parameters(ds::DynamicalSystem, ::Nothing) = merge(ds.dynamics, ds.securities)


# la llamo desde aca y no desde StochasticDiffEq? es mas liviano DiffEqBase
function DiffEqBase.SDEProblem(ds::DynamicalSystem{IIP}, tspan; u0=state(ds)) where {IIP}
    return SDEProblem{IIP}(
        SDEFunction(ds.f, ds.g),
        ds.g,
        u0,
        tspan,
        ds.params,
        noise=noise(ds),
        noise_rate_prototype=noise_rate_prototype(ds)
    )
end
