
# IDEA: implementar mejores dispatchs
#! SOLVE: es un temita el hecho de tener un dynamical system que resulte no ser DN pero que
#! tenga un modelo de short rate que es DN. Hay que solucionar eso

function drift!(
    μr::AbstractVector{<:Real}, r::Real, p::AffineParameters{OneFactor,true}, t::Real
)
    @unpack κ, θ = p
    μr[1] = κ(t) * (θ(t) - r)
    return nothing
end

drift(r::Real, p::AffineParameters{OneFactor,false}, t::Real) = p.κ(t) * (p.θ(t) - r)

function diffusion!(
  σr::AbstractVector{<:Real}, r::Real, p::AffineParameters{OneFactor,true}, t::Real
)
    @unpack Σ, α, β = p
    σr[1] = Σ(t) * sqrt(α(t) + β(t) * r)
    return nothing
end

diffusion(r::Real, p::AffineParameters{OneFactor,false}, t::Real) =
    p.Σ(t) * sqrt(p.α(t) + p.β(t) * r)

function drift!(
    μx::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    p::AffineParameters{MultiFactor,true},
    t::Real
)
    @unpack cache = p
    @unpack κ, θ, v1 = cache

    p.κ(κ, t)
    p.θ(θ, t)

    v1 .= θ .- x
    mul!(μx, κ, v1)

    return nothing
end

drift(x::AbstractVector{<:Real}, p::AffineParameters{MultiFactor,false}, t::Real) =
    p.κ(t) * (p.θ(t) - x)

# diagonal noise
function diffusion!(
    σx::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    p::AffineParameters{MultiFactor,true,D,true},
    t::Real
) where {D}
    @unpack cache = p
    @unpack Σ, α, β, v1 = cache

    p.Σ(Σ, t)
    p.α(α, t)
    p.β(β, t)

    # mul!(v1, β, x) .+= α(t)
    copyto!(v1, α)
    mul!(v1, β, x, true, true)
    map!(sqrt, v1, v1)
    σx .= Σ .* v1

    return nothing
end

# non-diagonal noise
function diffusion!(
    σx::AbstractMatrix{<:Real},
    x::AbstractVector{<:Real},
    p::AffineParameters{MultiFactor,true,D,false},
    t::Real
) where {D}

    @unpack cache = p
    @unpack Σ, α, β, v1 = cache

    p.Σ(Σ, t)
    p.α(α, t)
    p.β(β, t)

    # mul!(v1, β, x) .+= α(t)
    copyto!(v1, α)
    mul!(v1, β, x, true, true)
    map!(sqrt, v1, v1)
    mul!(σx, Σ, Diagonal(v1)) # permutedims is called and it allocates

    return nothing
end

# diagonal noise
function diffusion(
    x::AbstractVector{<:Real}, p::AffineParameters{MultiFactor,false,D,true}, t::Real
) where {D}
    @unpack Σ, α, β = p

    S = Diagonal(α(t) + β(t) * x)
    σ = Σ(t) * sqrt(S)

    # diag(::SMatrix) does not copy, otherwise we should use σ.diag
    return diag(σ)
end

# non-diagonal noise
function diffusion(
    x::AbstractVector{<:Real}, p::AffineParameters{MultiFactor,false,D,false}, t::Real
) where {D}
    @unpack Σ, α, β = p

    S = Diagonal(α(t) + β(t) * x)
    σ = Σ(t) * sqrt(S)

    # returns SMatrix since σ is non-diagonal (because Σ is not of Diagonal type)
    return σ
end

function drift!(
     μx::AbstractVector{<:Real}, x::Real, p::QuadraticParameters{OneFactor,true}, t::Real
)
    @unpack κ, θ = p
    μx[1] = κ(t) * (θ(t) - x)
    return nothing
end

drift(x::Real, p::QuadraticParameters{OneFactor,false}, t::Real) = p.κ(t) * (p.θ(t) - x)

function diffusion!(
    σx::AbstractVector{<:Real}, ::Real, p::QuadraticParameters{OneFactor,true}, t::Real
)
    @unpack σ = p
    σx[1] = σ(t)
    return nothing
end

diffusion(::Real, p::QuadraticParameters{OneFactor,false}, t::Real) = p.σ(t)

function drift!(
    μx::AbstractVector{<:Real},
    x::AbstractVector{<:Real},
    p::QuadraticParameters{MultiFactor,true},
    t::Real
)
    @unpack κ, θ = p
    @unpack v1 = p.cache

    v1 .= θ(t) .- x
    mul!(μx, κ(t), v1)

    return nothing
end

drift(x::AbstractVector{<:Real}, p::QuadraticParameters{MultiFactor,false}, t::Real) =
    p.κ(t) * (p.θ(t) - x)

# diagonal noise
function diffusion!(
    σx::AbstractVector{<:Real},
    ::AbstractVector{<:Real},
    p::QuadraticParameters{MultiFactor,true,D,true},
    t::Real
) where {D}
    @unpack σ = p
    copyto!(σx, σ(t).diag)
    return nothing
end

# non-diagonal noise
function diffusion!(
    σx::AbstractMatrix{<:Real},
    ::AbstractVector{<:Real},
    p::QuadraticParameters{MultiFactor,true,D,false},
    t::Real
) where {D}
    @unpack σ = p
    copyto!(σx, σ(t))
    return nothing
end

# diagonal noise
function diffusion(
    ::AbstractVector{<:Real}, p::QuadraticParameters{MultiFactor,false,D,true}, t::Real
) where {D}
    @unpack σ = p
    σt = σ(t) # es una matrix diagonal (de tipo Diagonal)
    return diag(σt) # diag(::SMatrix) does not copy, otherwise we should use σ.diag
end

# non-diagonal noise
function diffusion(
    ::AbstractVector{<:Real}, p::QuadraticParameters{MultiFactor,false,D,false}, t::Real
) where {D}
    @unpack σ = p
    σt = σ(t)
    return σt
end
