
struct HeathJarrowMortonModelInstantaneousForwardRate{H<:HeathJarrowMortonModelDynamics,S<:Security} <: InstantaneousForwardRate
    hjm::H
    f::S
end

(f::HeathJarrowMortonModelInstantaneousForwardRate)(t::Real, T::Real) = f.f(t, T)


"""
HeathJarrowMortonModelZeroCouponBond{L<:HeathJarrowMortonModelForwardRate} <: ZeroCouponBond

Zero coupon bond or discount bond obtained through the Heath Jarrow Morton Model.
"""
struct HeathJarrowMortonModelZeroCouponBond{T<:HeathJarrowMortonModelInstantaneousForwardRate} <: ZeroCouponBond
    f::T
end

"""
(P::HeathJarrowMortonModelZeroCouponBond)(t, T)

Computes and returns the zero coupon bond `P(t, T)` defined under the Heath Jarrow Morton Model. An
interpolation scheme is used in case `t` and/or `T` do not lie in the tenor structure of the
model.
"""
function (P::HeathJarrowMortonModelZeroCouponBond)(t::Real, T::Real)
    @unpack f = P
    @unpack Tenors, τ = get_parameters(P)

    if !(t in Tenors) || !(T in Tenors)
        throw(DomainError("`t` and `T` must lie in the tenor structure."))

    elseif t > T
        throw(DomainError("`t` must be ≤ `T`."))

    elseif 0 ≤ t < T
        return exp(-sum(f(i, t) * τ[i] for i = searchsortedfirst(Tenors, t):searchsortedfirst(Tenors, T)-1))

    elseif isequal(t, T)
        return one(Base.promote_eltype(1/t, 1/T))

    end
end

"""
HeathJarrowMortonModelForwardRate{f<:HeathJarrowMortonModel,P<:Security} <: SimpleForwardRate

Simply-compounded forward interest rate or simple forward rate obtained through the Heath Jarron Morton Model.
"""
struct HeathJarrowMortonModelForwardRate{P<:HeathJarrowMortonModelZeroCouponBond} <: SimpleForwardRate
    P::P
end


"""
(F::HeathJarrowMortonModelForwardRate)(t::Real, T::Real, S::Real)

Computes and returns the simple forward interest rate `F(t, T, S)`, with contiguous times
`T` and `S` over the Heath Jarrow Morton Model tenor structure.
"""
function (F::HeathJarrowMortonModelForwardRate)(t::Real, T::Real, S::Real)
    @unpack P = F
    @unpack Tenors = get_parameters(P)

    if !(t in Tenors) || !(T in Tenors) || !(S in Tenors)
        throw(DomainError("`t`, `T`, `S` must lie in the tenor structure."))

    elseif t > T
        throw(DomainError("`t` must be ≤ `T`."))

    elseif T >= S
        throw(DomainError("`T` must be < `S`."))

    elseif 0 ≤ t ≤ T < S
        return 1 / (S - T) * (P(t, T) / P(t, S) - 1)

    end
end

"""
(F::HeathJarrowMortonModelForwardRate)(T::Real, S::Real)

Computes and returns the simply-compounded spot interest rate, spot rate or Libor rate `F(T,
S) = F(T, T, S)`, with contiguous times `T` and `S` over the Heath Jarrow Morton Model tenor
structure.
"""
(F::HeathJarrowMortonModelForwardRate)(T::Real, S::Real) = F(T, T, S)


"""
HeathJarrowMortonModelMoneyMarketAccount{L<:HeathJarrowMortonModelForwardRate} <: DiscreteMoneyMarketAccount

Discrete money market account or bank account obtained through the Heath Jarrow Morton Model.
"""
struct HeathJarrowMortonModelMoneyMarketAccount{T<:HeathJarrowMortonModelInstantaneousForwardRate} <: DiscreteMoneyMarketAccount
    f::T
end

"""
(B::HeathJarrowMortonModelMoneyMarketAccount)(t::Real)

Computes and returns the money market account `B(t)` defined under the Heath Jarrow Morton Model.
An interpolation scheme is used in case `t` does not lie in the tenor structure of the
model.
"""
function (B::HeathJarrowMortonModelMoneyMarketAccount)(t::Real)
    @unpack f = B
    @unpack Tenors, τ = get_parameters(f)

    if !(t in Tenors)
        throw(DomainError("`t` must lie in the tenor structure."))

    elseif iszero(t)
        return one(1/t)

    elseif t > 0
       return exp(sum(f(i, t) * τ[i] for i = 1:searchsortedfirst(Tenors, t)-1))

    end

end

"""
HeathJarrowMortonModelDiscountFactor{B<:HeathJarrowMortonModelMoneyMarketAccount} <: DiscountFactor

Discount factor obtained through the Heath Jarrow Morton Model.
"""
struct HeathJarrowMortonModelDiscountFactor{T<:HeathJarrowMortonModelMoneyMarketAccount} <: DiscountFactor
    B::T
end

"""
(D::HeathJarrowMortonModelDiscountFactor)(t::Real, T::Real)

Computes and returns the discount factor `D(t, T)` defined under the Heath Jarrow Morton Model. An
interpolation scheme is used in case `t` and/or `T` do not lie in the tenor structure of the
model.
"""
function (D::HeathJarrowMortonModelDiscountFactor)(t::Real, T::Real)
    @unpack B = D
    @unpack Tenors = get_parameters(B)

    if !(t in Tenors) || !(T in Tenors)
        throw(DomainError("`t` and `T` must lie in the tenor structure."))

    elseif t > T
        throw(DomainError("`t` must be ≤ `T`."))

    elseif 0 ≤ t < T
        return B(t) / B(T)

    elseif isequal(t, T)
        return one(Base.promote_eltype(1/t, 1/T))
    end

end

struct HeathJarrowMortonModelSpotRate{T<:HeathJarrowMortonModelInstantaneousForwardRate} <: SpotRate
    f::T
end

function (r::HeathJarrowMortonModelSpotRate)(t::Real)
    @unpack f = r
    @unpack Tenors = get_parameters(f)

    if !(t in Tenors)
        throw(DomainError("`t` must lie in the tenor structure."))

    elseif t > 0
        return f(searchsorted(Tenors, t)[1], t)

    elseif iszero(t)
        return f(1)
    end

end


function FixedIncomeSecurities(hjm::HJM, f::Security) where {HJM<:HeathJarrowMortonModelDynamics}
    f = HeathJarrowMortonModelInstantaneousForwardRate(hjm, f)
    P = HeathJarrowMortonModelZeroCouponBond(f)
    B = HeathJarrowMortonModelMoneyMarketAccount(f)
    D = HeathJarrowMortonModelDiscountFactor(B)
    r = HeathJarrowMortonModelSpotRate(f)
    F = HeathJarrowMortonModelForwardRate(P)

    return FixedIncomeSecurities{HJM}(r, B, D, P, F, f)
end

get_parameters(f::HeathJarrowMortonModelInstantaneousForwardRate) = get_parameters(f.hjm)
get_parameters(P::HeathJarrowMortonModelZeroCouponBond) = get_parameters(P.f.hjm)
get_parameters(B::HeathJarrowMortonModelMoneyMarketAccount) = get_parameters(B.f.hjm)


