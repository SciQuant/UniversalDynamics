struct Schlogl{A} <: LiborMarketModelInterpolationMethod
    Schlogl(α::Bool = true) = new{α}()
end

# TODO: interpolation de las L(t, T, S), deberia ser posible

# esto es con α(t) = 1
"""
    interpolate(P::LiborMarketModelZeroCouponBond, ::Schlogl{<:Nothing}, t::Real, T::Real)

Computes the zero coupon bond price `P(t, T)` when `t` and/or `T` do not lie in the tenor
structure of the Libor Market Model using the Schlogl interpolation method assuming zero
volatility for short dated bonds.
"""
function interpolate(P::LiborMarketModelZeroCouponBond, ::Type{Schlogl{false}}, t::Real, T::Real)
    @unpack L = P
    @unpack Tenors, τ = L.params

    # FIXME: estas busquedas de alguna forma ya las hicimos en rt y rT en la funcion
    # anterior, evitarlas! basicamente es pensar la relacion entre rt y qt si llegue hasta
    # esta funcion y lo mismo para qT y rT. Debugeando salen faciles, de hecho creo que ya
    # lo hice.
    qt = q(Tenors, t)
    qT = q(Tenors, T)

    # short-dated bond
    Pttq = 1 / (1 + (Tenors[qt] - t) * L(qt - 1, Tenors[qt - 1]))

    E_1_PTTq = 1 + (Tenors[qT] - T) * L(qT - 1, t)

    range = qt:qT-1
    if isempty(range)
        return Pttq * E_1_PTTq
    else
        PtqTq = prod((1 / (1 + τ[i] * L(i, t))) for i in range)
        return Pttq * PtqTq * E_1_PTTq
    end
end

α(Tenors::Vector{<:Real}, t::Real) =  α(Tenors, t, q(Tenors, t))
α(Tenors::Vector{<:Real}, t::Real, qt::Int64) = (Tenors[qt] - t) / (Tenors[qt] - Tenors[qt - 1])

"""
    interpolate(P::LiborMarketModelZeroCouponBond, imethod::Schlogl{F}, t::Real, T::Real) where {F<:Function}

Computes the zero coupon bond price `P(t, T)` when `t` and/or `T` do not lie in the tenor
structure of the Libor Market Model using the Schlogl interpolation method relaxing the zero
volatility for short dated bonds... describir α(t).
"""
function interpolate(P::LiborMarketModelZeroCouponBond, ::Type{Schlogl{true}}, t::Real, T::Real)
    @unpack L = P
    @unpack Tenors, τ, σ = parameters(P)

    # notar que si T cae en el ultimo intervalo de la Tenor, llamo una L que no esta
    # simulada, lo que implica un error. Lo mismo ocurre si t tambien cae en el ultimo
    # intervalo (y en consecuencia, T tambien). Todos estos casos requieren la simulacion de
    # una Libor mas, que expire en la ultima fecha de la tenor y madure mas adelante. Sin
    # embargo, no le encuentro sentido a esto, por lo que no lo tendre en cuenta y
    # considerare que T < Tenors[end-1].

    qt = q(Tenors, t)
    qT = q(Tenors, T)

    αt = α(Tenors, t, qt)
    αT = α(Tenors, T, qT)

    Pttq = 1 / (1 + (Tenors[qt] - t) * (αt * L(qt - 1, Tenors[qt - 1]) + (1 - αt) * L(qt, t)))

    # correction factor
    LqTt = L(qT, t)

    #! now it should be (σ(s)[qt])^2 and use cache for IIP, ademas esto es para DN porque es
    #! 1 componente no nula. Cuando σ retorna un vector, ver la integral de Proposition 2.2
    #! paper copado. Ahi se ve que λ(t, Tᵢ) retorna un vector fila y tenemos λ ρ λᵀ. Habria
    #! que ver que es ρ, si la matrix o una componente. Pero casi seguro es la matrix. Si,
    #! es la matrix.
    int = integral(s -> σ(qT, s)^2, t, T)

    Δ = 1 + ((Tenors[qT + 1] - Tenors[qT]) * LqTt * (exp(int) - 1)) / (1 + (Tenors[qT + 1] - Tenors[qT]) * LqTt)

    E_1_PTTq = 1 + (Tenors[qT] - T) * (αT * L(qT - 1, t) + (1 - αT) * LqTt * Δ)

    range = qt:qT-1
    if isempty(range)
        return Pttq * E_1_PTTq
    else
        PtqTq = prod((1 / (1 + τ[i] * L(i, t))) for i in range)
        return Pttq * PtqTq * E_1_PTTq
    end
end

# CHECK: esta funcion no la pense mucho, simplemente use el P(t, T_{q(t)}) que calculo en
# interpolate(P::LiborMarketModelZeroCouponBond, ...)
"""
    interpolate(B::LiborMarketModelMoneyMarketAccount, ::Schlogl{<:Nothing}, t::Real)

Computes the money market account value `B(t)` when `t` does not lie in the tenor structure
of the Libor Market Model using the Schlogl interpolation method assuming zero volatility
for short dated bonds.
"""
function interpolate(B::LiborMarketModelMoneyMarketAccount, ::Type{Schlogl{false}}, t::Real)
  @unpack L = B
  @unpack Tenors, τ = L.params

  # NOTE: este metodo sirve aun cuando t ∈ Tenors. La unica diferencia con lo de la funcion
  # madre es que me evito hacer unas cuentas demas. Ver si unimos todo.

  # FIXME: esta busqueda de alguna forma ya la hice en r en la funcion anterior, evitarla!
  # basicamente es pensar la relacion entre r y qt si llegue hasta esta funcion.
  qt = q(Tenors, t)

  # short-dated bond
  Pttq = 1 / (1 + (Tenors[qt] - t) * L(qt - 1, Tenors[qt - 1]))

  # este range nunca puede estar empty
  Btq = prod((1 + τ[i] * L(i, t)) for i in 1:qt-1)

  return Pttq * Btq
end