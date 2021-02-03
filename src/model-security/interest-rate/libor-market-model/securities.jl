
"""
    LiborMarketModelForwardRate{L<:LiborMarketModel,P<:Security} <: SimpleForwardRate

Simply-compounded forward interest rate or simple forward rate obtained through the Libor
Market Model.
"""
struct LiborMarketModelForwardRate{L<:LiborMarketModelDynamics,P<:Security} <: SimpleForwardRate
    lmm::L
    Ln::P
end

"""
    (L::LiborMarketModelForwardRate)(t::Real)

Computes and returns all the simple forward interest rates at time `t` described by the
implemented Libor Market Model such that `L_n(t) = L(t, T_n, T_{n+1})` is computed for each
`n` and with `T` the Libor Market Model tenor structure.
"""
(L::LiborMarketModelForwardRate)(t::Real) = L.Ln(t)

"""
    (L::LiborMarketModelForwardRate)(i::Int64, t::Real)

Computes and returns the `i`-th simple forward interest rates at time `t` described by the
implemented Libor Market Model such that `L_n(t) = L(t, T_n, T_{n+1})` is computed
for `n = i` and with `T` the Libor Market Model tenor structure.
"""
(L::LiborMarketModelForwardRate)(i::Int64, t::Real) = L.Ln(i, t)

"""
    (L::LiborMarketModelForwardRate)(idxs::Union{Vector{Int64},StepRange}, t::Real)

Computes and returns some simple forward interest rates at time `t` described by the
implemented Libor Market Model such that `L_n(t) = L(t, T_n, T_{n+1})` is computed for `n in
idxs` and with `T` the Libor Market Model tenor structure.
"""
(L::LiborMarketModelForwardRate)(idxs::Union{Vector{Int64},StepRange}, t::Real) = L.Ln(idxs, t)

"""
    (L::LiborMarketModelForwardRate)(t::Real, T::Real, S::Real)

Computes and returns the simple forward interest rate `L(t, T, S)`, with contiguous times
`T` and `S` over the Libor Market Model tenor structure. In case `T` and `S` are not
contiguous... (we should investigate this and allow for a interpolation scheme).
"""
function (L::LiborMarketModelForwardRate)(t::Real, T::Real, S::Real)
  @unpack Tenors = parameters(L)

    if 0 ≤ t ≤ T < S
        rT = searchsorted(Tenors, T)
        rS = searchsorted(Tenors, S)

        if !isempty(rT) && !isempty(rS)
            i = getindex(rT)
            j = getindex(rS)

            if (i + 1) != j
                # TODO: aca podria existir la interpolacion que discuto abajo?
                throw(DomainError("continuous..."))
            else
                return L(i, t)
            end
        else
            # creo q con el paper copado podria implementar un metodo de interpolacion, ver
            # fig 2.3 en teoria incluso si tengo los P's interpolados, las L's deberian
            # poder interpolarse
            # TODO: return interpolate(L, t, T)
            throw(DomainError("Forward interest rate interpolation not yet implemented."))
        end
    elseif isequal(T, S)
        throw(DomainError("error 1..."))
    elseif S < T
        throw(DomainError("error 2..."))
    else
        throw(DomainError("error 3..."))
    end
end

"""
    (L::LiborMarketModelForwardRate)(T::Real, S::Real)

Computes and returns the simply-compounded spot interest rate, spot rate or Libor rate `L(T,
S) = L(T, T, S)`, with contiguous times `T` and `S` over the Libor Market Model tenor
structure. In case `T` and `S` are not contiguous... (see above).
"""
(L::LiborMarketModelForwardRate)(T::Real, S::Real) = L(T, T, S)

"""
    LiborMarketModelZeroCouponBond{L<:LiborMarketModelForwardRate} <: ZeroCouponBond

Zero coupon bond or discount bond obtained through the Libor Market Model.
"""
struct LiborMarketModelZeroCouponBond{T<:LiborMarketModelForwardRate} <: ZeroCouponBond
    L::T
end

"""
    (P::LiborMarketModelZeroCouponBond)(t, T)

Computes and returns the zero coupon bond `P(t, T)` defined under the Libor Market Model. An
interpolation scheme is used in case `t` and/or `T` do not lie in the tenor structure of the
model.
"""
function (P::LiborMarketModelZeroCouponBond)(t::Real, T::Real)
    @unpack L = P
    @unpack Tenors, τ, imethod = parameters(P)

    if 0 ≤ t < T
        rt = searchsorted(Tenors, t)
        rT = searchsorted(Tenors, T)

        # Check if t and T belong to Tenors assuming they are not repeated. Since we
        # constructed the Tenors object, we could relax the last check and just check if the
        # ranges are not empty.
        if !isempty(rt) && !isempty(rT)
            # usamos la interpretacion q(Tn) = n, i.e. left continuous, para este calculo
            n = getindex(rt)
            N = getindex(rT)
            return prod((1 / (1 + τ[i] * L(i, t))) for i in n:N-1)
        else
            return interpolate(P, imethod, t, T)
        end
    elseif isequal(t, T)
        return one(Base.promote_eltype(1/t, 1/T))
    else
        throw(DomainError("`t` must be ≤ `T` when computing a Zero Coupon Bond P(t, T)."))
    end
end

"""
    LiborMarketModelMoneyMarketAccount{L<:LiborMarketModelForwardRate} <: DiscreteMoneyMarketAccount

Discrete money market account or bank account obtained through the Libor Market Model.
"""
struct LiborMarketModelMoneyMarketAccount{T<:LiborMarketModelForwardRate} <: DiscreteMoneyMarketAccount
    L::T
end

"""
    (B::LiborMarketModelMoneyMarketAccount)(t::Real)

Computes and returns the money market account `B(t)` defined under the Libor Market Model.
An interpolation scheme is used in case `t` does not lie in the tenor structure of the
model.
"""
function (B::LiborMarketModelMoneyMarketAccount)(t::Real)
    @unpack L = B
    @unpack Tenors, τ, imethod = parameters(B)

    if t > 0
        r = searchsorted(Tenors, t)
        if !isempty(r)
            # usamos la interpretacion q(Tn) = n, i.e. left continuous, para este calculo
            n = getindex(r)
            return prod((1 + τ[i] * L(i, t)) for i in 1:n-1)
        else
            return interpolate(B, imethod, t)
        end
    elseif iszero(t)
        return one(1/t)
    else
        throw(DomainError("error?"))
    end
end

"""
    LiborMarketModelDiscountFactor{B<:LiborMarketModelMoneyMarketAccount} <: DiscountFactor

Discount factor obtained through the Libor Market Model.
"""
struct LiborMarketModelDiscountFactor{T<:LiborMarketModelMoneyMarketAccount} <: DiscountFactor
    B::T
end

"""
    (D::LiborMarketModelDiscountFactor)(t::Real, T::Real)

Computes and returns the discount factor `D(t, T)` defined under the Libor Market Model. An
interpolation scheme is used in case `t` and/or `T` do not lie in the tenor structure of the
model.
"""
function (D::LiborMarketModelDiscountFactor)(t::Real, T::Real)
    @unpack B = D
    if 0 ≤ t < T
        return B(t) / B(T)
    elseif isequal(t, T)
        return one(Base.promote_eltype(1/t, 1/T))
    else
        throw(DomainError("`t` must be ≤ `T` when computing a Discount Factor D(t, T)."))
    end
end

struct LiborMarketModelSpotRate <: SpotRate end

(r::LiborMarketModelSpotRate)(args...) =
    error("the spot rate is not defined in the Libor Market Model.")

struct LiborMarketModelInstantaneousForwardRate <: InstantaneousForwardRate end

(r::LiborMarketModelInstantaneousForwardRate)(args...) =
    error("the instantaneous forward rate is not defined in the Libor Market Model.")

function FixedIncomeSecurities(lmm::LMM, Ln::Security) where {LMM<:LiborMarketModelDynamics}
    L = LiborMarketModelForwardRate(lmm, Ln)
    r = LiborMarketModelSpotRate()
    B = LiborMarketModelMoneyMarketAccount(L)
    D = LiborMarketModelDiscountFactor(B)
    P = LiborMarketModelZeroCouponBond(L)
    f = LiborMarketModelInstantaneousForwardRate()
    return FixedIncomeSecurities{LMM}(r, B, D, P, L, f)
end

parameters(L::LiborMarketModelForwardRate) = parameters(L.lmm)
parameters(P::LiborMarketModelZeroCouponBond) = parameters(P.L.lmm)
parameters(B::LiborMarketModelMoneyMarketAccount) = parameters(B.L.lmm)