
"""
    BasicFixedIncomeSecurity

Supertype for the basic fixed income market securities which are fundamental for Interest
Rate modeling.
"""
abstract type BasicFixedIncomeSecurity <: ModelSecurity end

"""
    SpotRate <: BasicFixedIncomeSecurity

Supertype for spot rates, also known as short rates.
"""
abstract type SpotRate <: BasicFixedIncomeSecurity end

"""
    MoneyMarketAccount <: BasicFixedIncomeSecurity

Supertype for money market accounts, also known as bank accounts.
"""
abstract type MoneyMarketAccount <: BasicFixedIncomeSecurity end

"""
    ContinuousMoneyMarketAccount <: MoneyMarketAccount

Supertype for continuous money market accounts ``\beta(t, T)``.
"""
abstract type ContinuousMoneyMarketAccount <: MoneyMarketAccount end

"""
    DiscreteMoneyMarketAccount <: MoneyMarketAccount

Supertype for discrete money market accounts ``B(t, T)``.
"""
abstract type DiscreteMoneyMarketAccount <: MoneyMarketAccount end

"""
    DiscountFactor <: BasicFixedIncomeSecurity

Supertype for discount factors ``D(t, T)``.
"""
abstract type DiscountFactor <: BasicFixedIncomeSecurity end

"""
    ZeroCouponBond <: BasicFixedIncomeSecurity

Supertype for zero coupon bonds, also known as discount bonds, ``P(t, T)``.
"""
abstract type ZeroCouponBond <: BasicFixedIncomeSecurity end

"""
    SimpleForwardRate <: BasicFixedIncomeSecurity

Supertype for the simply-compounded forward interest rate or simple forward rate ``L(t, T,
S)``.
"""
abstract type SimpleForwardRate <: BasicFixedIncomeSecurity end

"""
    InstantaneousForwardRate <: BasicFixedIncomeSecurity

Supertype for the instantaneous forward rate ``f(t, T)``.
"""
abstract type InstantaneousForwardRate <: BasicFixedIncomeSecurity end

"""
    FixedIncomeSecurities

Defines an insterest rate object with basic fixed income securities given by an specific
interest rate model, such as the [`ShortRateModel`](@ref) or the [`LiborMarketModel`](@ref).
"""
struct FixedIncomeSecurities{
    MD<:InterestRateModelDynamics,
    T1<:SpotRate,
    T2<:MoneyMarketAccount,
    T3<:DiscountFactor,
    T4<:ZeroCouponBond,
    T5<:SimpleForwardRate,
    T6<:InstantaneousForwardRate
}
    r::T1
    B::T2
    D::T3
    P::T4
    L::T5
    f::T6
end

FixedIncomeSecurities{MD}(r::T1, B::T2, D::T3, P::T4, L::T5, f::T6) where {MD,T1,T2,T3,T4,T5,T6} =
    FixedIncomeSecurities{MD,T1,T2,T3,T4,T5,T6}(r, B, D, P, L, f)

include("interest-rate/short_rate_model.jl")
# include("interest-rate/libor_market_model.jl")
