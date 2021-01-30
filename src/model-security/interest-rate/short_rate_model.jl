include("riccati.jl")


struct ShortRateModelSpotRate{S<:ShortRateModelDynamics,X,R} <: SpotRate
    srm::S
    x::X
    r::R
end

function ShortRateModelSpotRate(srm::OneFactorAffineModelDynamics, r::Security)
    @unpack ξ₀, ξ₁ = parameters(srm)
    x(t::Real) = (r(t) - ξ₀(t)) / ξ₁(t)
    return ShortRateModelSpotRate(srm, x, r)
end

function ShortRateModelSpotRate(srm::MultiFactorAffineModelDynamics{false}, x::Security)
    @unpack ξ₀, ξ₁ = parameters(srm)
    r(t::Real) = ξ₀(t) + dot(ξ₁(t), x(t))
    return ShortRateModelSpotRate(srm, x, r)
end

function ShortRateModelSpotRate(srm::MultiFactorAffineModelDynamics{true}, x::Security)
    p = parameters(srm)
    @unpack cache = p
    @unpack ξ₁ = cache
    function r(t::Real)
        ξ₀ = p.ξ₀(t)
        p.ξ₁(ξ₁, t)
        return ξ₀ + dot(ξ₁, x(t))
    end
    return ShortRateModelSpotRate(srm, x, r)
end

function ShortRateModelSpotRate(srm::OneFactorQuadraticModelDynamics, x::Security)
    @unpack ξ₀, ξ₁, ξ₂ = parameters(srm)

    r = function (t::Real)
        xt = x(t)
        return ξ₀(t) + ξ₁(t) * xt + ξ₂(t) * xt^2
        # return @muladd ξ₀(t) + ξ₁(t) * xt + ξ₂(t) * xt^2
    end

    return ShortRateModelSpotRate(srm, x, r)
end

function ShortRateModelSpotRate(srm::MultiFactorQuadraticModelDynamics, x::Security)
    @unpack ξ₀, ξ₁, ξ₂ = parameters(srm)

    r = function (t::Real)
        xt = x(t)
        return ξ₀(t) + dot(ξ₁(t), xt) + dot(xt, ξ₂(t), xt)
    end

    return ShortRateModelSpotRate(srm, x, r)
end

(r::ShortRateModelSpotRate)(t::Real)= r.r(t)

struct ShortRateModelMoneyMarketAccount{T} <: ContinuousMoneyMarketAccount
    B::T
end

(B::ShortRateModelMoneyMarketAccount)(t::Real) = B.B(t)

struct ShortRateModelDiscountFactor{T<:ShortRateModelMoneyMarketAccount} <: DiscountFactor
    B::T
end

# esta funcion es igual a LMMDiscountFactor asi que podriamos hacerla para D::DiscountFactor
function (D::ShortRateModelDiscountFactor)(t::Real, T::Real)
    @unpack B = D
    if 0 ≤ t < T
        return B(t) / B(T)
    elseif isequal(t, T)
        return one(Base.promote_eltype(1/t, 1/T))
    else
        throw(DomainError("`t` must be ≤ `T` for Discount Factor D(t, T)."))
    end
end

struct ShortRateModelZeroCouponBond{R<:ShortRateModelSpotRate} <: ZeroCouponBond
    r::R
end

shortratemodel(P::ShortRateModelZeroCouponBond) = P.r.srm

function (P::ShortRateModelZeroCouponBond)(t::Real, T::Real, xt::Union{Real,AbstractVector{<:Real}})
    if 0 ≤ t < T
        return zerocouponbond(shortratemodel(P), P, t, T, xt) # dispatch by model
    elseif isequal(t, T)
        return one(eltype(P.r.srm)) # one(Base.promote_eltype(1/t, 1/T))
    else
        throw(DomainError("`t` must be ≤ `T` for Zero Coupon Bond P(t, T)."))
    end
end

# TODO: allocates x(t), use x(out, t) if Array
(P::ShortRateModelZeroCouponBond)(t::Real, T::Real) = P(t, T, P.r.x(t))

# TODO: las allocations de las funciones zerocouponbond hay que estudiarlas para cuando nos
# movemos con SArrays o Arrays, ya que las estudie muy poco y seguramente puede mejorarse el
# asunto. Tengo que usar la variable rout del cache en el modelo.

function zerocouponbond(
    ::OneFactorAffineModelDynamics, P::ShortRateModelZeroCouponBond, t::Real, T::Real, xt::Real
)
    u = solve_riccati(P, t, T)
    a, b = u
    # aca falta multiplicar por \frac{\exp \left( \int_0^T \xi_0(s) ds \right)}{\exp \left( \int_0^t \xi_0(s) ds \right}.
    # Podemos pedir un parametro que sea P(s) = \exp \left( \int_0^s \xi_0(s) \right) ya que
    # comunmente ξ₀(t) = fᴹ(0, t), por lo que P(s) = Pᴹ(0, s) y evitamos hacer las
    # integrales. numericamente. Otra idea que habria que experimentar es ver si podemos
    # subir este termino al exponente y que quede en funcion de rt? O sea, que aparezca algo
    # como x(t) - ξ₀(t)?
    return exp(a - b * xt)
    # @muladd return exp(a - B * xt)
end

function zerocouponbond(
    ::MultiFactorAffineModelDynamics,
    P::ShortRateModelZeroCouponBond,
    t::Real,
    T::Real,
    xt::AbstractVector{<:Real}
)
    # TODO: esta funcion hay que analizarla si `u` es un SVector
    u = solve_riccati(P, t, T) # aca seguro alloca si recibo un vect? ver para IIP
    a = u[1]
    b = view(u, 2:lastindex(u))
    return exp(a - dot(b, xt))
end

function zerocouponbond(
    ::OneFactorQuadraticModelDynamics, P::ShortRateModelZeroCouponBond, t::Real, T::Real, xt::Real
)
    u = solve_riccati(P, t, T)
    a, b, c = u
    return exp(-a - b * xt - c * xt^2)
    # @muladd return exp(-a - b * xt - c * xt^2)
end

function zerocouponbond(
    srm::MultiFactorQuadraticModelDynamics,
    P::ShortRateModelZeroCouponBond,
    t::Real,
    T::Real,
    xt::AbstractVector{<:Real}
)
    # TODO: esta funcion hay que analizarla si `u` es un SVector, ya que luego los view y
    # reshape puede que anden mal. Quizas hay que usar el MultiFactorQuadratic{FM,D} y tomar
    # de aca el D y hacer Size(D, D), todo si es OOP. Lo mismo hay que ver para los otros
    # casos
    u = solve_riccati(P, t, T)
    N = dimension(srm)
    a = u[1]
    b = view(u, 2:N+1)
    c⃗ = view(u, N+2:lastindex(u))
    c = reshape(c⃗, (N, N))
    return exp(-a - dot(b, xt) - dot(xt, c, xt))
end

function solve_riccati(P::ShortRateModelZeroCouponBond, t::Real, T::Real)
    # solve Riccati only if needed
    srm = shortratemodel(P)
    sol = remake_and_solve(srm.prob, T)

    # for now, both allocate (see #1270 in OrdinaryDiffEq)
    if isinplace(srm)
        u = srm.cache.rout
        sol(u, t)
    else
        u = sol(t)
    end
    return u # mmm
end

# Las "callables" aun no son thread safe, incluso usando ThreadSafeDict al parecer. Por lo
# tanto, habra que ver si tenemos que trabajar con Memoize y LRU caches. Aca estamos usando
# como cache el default. Quiero destacar que esta tarda casi lo mismo con Memoize que con
# Memoization y es el mismo codigo, pero Memoization tenia una allocation
Memoize.@memoize function remake_and_solve(prob, T)
    println("Running")
    newprob = OrdinaryDiffEq.remake(prob, tspan = (T, zero(eltype(prob.u0))))
    sol = OrdinaryDiffEq.solve(newprob, Tsit5())
    return sol
end

# esta funciona perfectamente pero para el caso de estudio era mucho mas lenta que la
# anterior y ademas alloca. Para Memoization no pude poner todo el dict con types porque no
# lo tomaba bien. Cuando se hace el @expand, se ve que es raro y la funcion anterior si
# hace lo que uno esperaria hacer!
# Memoize.@memoize Dict{Tuple{ODEProblem,Real},ODESolution} function remake_and_solve(prob, T)
#     println("Running")
#     newprob = remake(prob, tspan = (T, zero(eltype(prob.u0))))
#     sol = OrdinaryDiffEq.solve(newprob, Tsit5())
#     return sol
# end

struct ShortRateModelForwardRate{T<:ShortRateModelZeroCouponBond} <: SimpleForwardRate
    P::T
end

function (L::ShortRateModelForwardRate)(t::Real, T::Real, S::Real)
    P = L.P
    return 1 / (S - T) * (P(t, T) / P(t, S) - 1)
end

(L::ShortRateModelForwardRate)(T::Real, S::Real) = L(T, T, S)

struct ShortRateModelInstantaneousForwardRate{
    T<:ShortRateModelZeroCouponBond
} <: InstantaneousForwardRate
    P::T
end

# ver si cambia la cosa si paso el closure a la estructura de InstantaneousForwardRate
(f::ShortRateModelInstantaneousForwardRate)(t::Real, T::Real) = ForwardDiff.derivative(s -> -log(f.P(t, s)), T)

function FixedIncomeSecurities(srmd::SRMD, x::Security, β::Security) where {SRMD<:ShortRateModelDynamics}
    r = ShortRateModelSpotRate(srmd, x)
    B = ShortRateModelMoneyMarketAccount(β)
    D = ShortRateModelDiscountFactor(B)
    P = ShortRateModelZeroCouponBond(r)
    L = ShortRateModelForwardRate(P)
    f = ShortRateModelInstantaneousForwardRate(P)
    return FixedIncomeSecurities{SRMD}(r, B, D, P, L, f)
end
