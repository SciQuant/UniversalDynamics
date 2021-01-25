
"""
    AbstractMeasure

Supertype for all measures.
"""
abstract type AbstractMeasure end

@doc raw"""
    RiskNeutral <: AbstractMeasure

Also known as ``Q``, its numeraire is the continuously compounded money market account
``\beta(t)``.
"""
struct RiskNeutral <: AbstractMeasure end

@doc raw"""
    Spot <: AbstractMeasure

Uses the discrete-time equivalent of the continuously compounded money market account as the
numeraire. It is often convenient when working with a multitude of forward rates on a tenor
structure.
"""
struct Spot <: AbstractMeasure end

# ojo, esta se puede mal interpretar como Tn, que es una measure variable, y no es lo que
# queremos representar aqui.
@doc raw"""
    TForward{I} <: AbstractMeasure

Also known as ``Q^T``, uses a ``T`` maturity zero coupon bond ``P(t, T)`` as the numeraire.
When working with a tenor structure, `I` identifies ``P(t, T_I)`` as the numeraire, with
``T_I`` the `I`-th maturity in the tenor structure.
"""
struct TForward{N} <: AbstractMeasure end

@doc raw"""
    Terminal <: AbstractMeasure

Uses a ``T_N`` maturity zero coupon bond ``P(t, T_N)``, with ``T_N`` the last maturity in
the tenor structure, as numeraire.
"""
struct Terminal <: AbstractMeasure end
