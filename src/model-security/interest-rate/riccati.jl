
riccati_dimension(::Type{<:OneFactorAffineModelDynamics}) = 2
riccati_dimension(::Type{<:MultiFactorAffineModelDynamics{IIP,D}}) where {IIP,D} = 1 + D
riccati_dimension(::Type{<:OneFactorQuadraticModelDynamics}) = 3
riccati_dimension(::Type{<:MultiFactorQuadraticModelDynamics{IIP,D}}) where {IIP,D} = 1 + D * (1 + D)

function riccati_problem(
    T::Type{<:ShortRateModelDynamics{FM,IIP}}, ::Type{S}, p
) where {FM,IIP,S}
    N = riccati_dimension(T)
    uT = IIP ? zeros(S, N) : SVector{N}(zeros(S, N))
    prob =  ODEProblem{IIP}(riccati, uT, (one(S), zero(S)), p)
    return prob
end

function riccati(du, u, p::AffineParameters{OneFactor,true}, t)
    κ, θ, Σ, α, β, ξ₀, ξ₁ = p(t)

    ξ₀′ = ForwardDiff.derivative(p.ξ₀, t)
    ξ₁′ = ForwardDiff.derivative(p.ξ₁, t)

    @inbounds begin
        B = u[2]

        # use @muladd here? check its expansion
        # @muladd begin
        du[1] = (κ * θ - κ * ξ₀ - ξ₀′) * B / ξ₁ - (α + β * ξ₀) / 2 * (Σ * B / ξ₁)^2
        du[2] = (κ + ξ₁′ / ξ₁) * B + β / (2 * ξ₁) * (Σ * B)^2 - ξ₁
    end

  return nothing
end

function riccati(u, p::AffineParameters{OneFactor,false}, t)
    κ, θ, Σ, α, β, ξ₀, ξ₁ = p(t)

    ξ₀′ = ForwardDiff.derivative(p.ξ₀, t)
    ξ₁′ = ForwardDiff.derivative(p.ξ₀, t)

    @inbounds  B = u[2]

    # use @muladd here? check its expansion
    # @muladd begin
    dA = (κ * θ - κ * ξ₀ - ξ₀′) * B / ξ₁ - (α + β * ξ₀) / 2 * (Σ * B / ξ₁)^2
    dB = (κ + ξ₁′ / ξ₁) * B + β / (2 * ξ₁) * (Σ * B)^2 - ξ₁

    return SVector{2}(dA, dB)
end

function riccati(du, u, p::AffineParameters{MultiFactor,true}, t)
    κ, θ, Σ, α, β, ξ₀, ξ₁ = p(t)
    @unpack v1, v2, v3 = p.cache

    @inbounds begin
        B = @view u[2:end]

        # TODO: reducir caches
        mul!(v1, transpose(κ), B)
        mul!(v2, transpose(Σ), B)
        mul!(v3, Diagonal(v2), v2)

        du[1] = dot(θ, v1) - dot(α, v3) / 2 + ξ₀
        du[2:end] .= v1 .+ mul!(v2, transpose(β), v3, 1/2, 0) .- ξ₁
    end

    return nothing
end

# TODO: in place version with `DN = true`, i.e. Σ returns a vector instead of a matrix.

function riccati(u, p::AffineParameters{MultiFactor,false}, t)
    κ, θ, Σ, α, β, ξ₀, ξ₁ = p(t)

    B = @view u[2:end]

    κᵀB = transpose(κ) * B
    ΣᵀB = transpose(Σ) * B
    diagΣᵀB_ΣᵀB = Diagonal(ΣᵀB) * ΣᵀB

    dA = dot(θ, κᵀB) - dot(α, diagΣᵀB_ΣᵀB) / 2 + ξ₀
    dB = κᵀB + transpose(β) * diagΣᵀB_ΣᵀB / 2 - ξ₁

    return vcat(SVector(dA), dB)
end

function riccati(u, p::AffineParameters{MultiFactor,false,D,true}, t) where {D}
    κ, θ, Σ, α, β, ξ₀, ξ₁ = p(t)

    B = @view u[2:end]

    κᵀB = transpose(κ) * B
    ΣᵀB = Σ .* B
    diagΣᵀB_ΣᵀB = ΣᵀB .* ΣᵀB

    dA = dot(θ, κᵀB) - dot(α, diagΣᵀB_ΣᵀB) / 2 + ξ₀
    dB = κᵀB + transpose(β) * diagΣᵀB_ΣᵀB / 2 - ξ₁

    return vcat(SVector(dA), dB)
end

function riccati(du, u, p::QuadraticParameters{OneFactor,true}, t)
    κ, θ, σ, ξ₀, ξ₁, ξ₂ = p(t)

    B = u[2]
    C = u[3]

    # use @muladd here? check its expansion
    # @muladd begin
    du[1] = - σ^2 * C - κ * θ * B + (σ * B)^2 / 2 - ξ₀
    du[2] = κ * B + 2 * σ^2 * B * C - 2 * κ * θ * C - ξ₁
    du[3] = 2 * κ * C + 2 * (σ * C)^2 - ξ₂

    return nothing
end

function riccati(u, p::QuadraticParameters{OneFactor,false}, t)
    κ, θ, σ, ξ₀, ξ₁, ξ₂ = p(t)

    B = u[2]
    C = u[3]

    # use @muladd here? check its expansion
    # @muladd begin
    dA = - σ^2 * C - κ * θ * B + (σ * B)^2 / 2 - ξ₀
    dB = κ * B + 2 * σ^2 * B * C - 2 * κ * θ * C - ξ₁
    dC = 2 * κ * C + 2 * (σ * C)^2 - ξ₂

    return SVector{3}(dA, dB, dC)
end

function riccati(du, u, p::QuadraticParameters{MultiFactor,true,D}, t) where {D}
    κ, θ, σ, ξ₀, ξ₁, ξ₂ = p(t)

    @unpack v1, v2, m1, m2 = p.cache

    Brange = 2:D+1
    Crange = D+2:lastindex(u)

    B = @view u[Brange]
    C⃗ = @view u[Crange]
    C = reshape(C⃗, (D, D))

    κᵀ = transpose(κ)
    mul!(v1, κ, θ)
    mul!(m1, σ, transpose(σ)) # allocates for N = 2 or N = 3, issue #37873 (JuliaLang/julia)
    mul!(v2, m1, B)
    mul!(m2, m1, C)

    du[1] = -tr(m2) - dot(B, v1) + dot(B, v2) / 2 - ξ₀

    dB = view(du, Brange)
    mul!(dB, κᵀ, B,  1, 0)
    mul!(dB, C, v2,  2, 1)
    mul!(dB, C, v1, -2, 1)
    dB .-= ξ₁

    mul!(m1, transpose(C), m2, 2, 0)
    mul!(m1, κᵀ, C, true, true)
    mul!(m1, C, κ, true, true)
    m1 .-= ξ₂

    copyto!(view(du, Crange), m1)

    return nothing
end

function riccati(u, p::QuadraticParameters{MultiFactor,false,D}, t) where {D}
    κ, θ, σ, ξ₀, ξ₁, ξ₂ = p(t)

    Brange = 2:D+1
    Crange = D+2:lastindex(u)

    B = @view u[Brange]
    C⃗ = @view u[Crange]
    C = reshape(C⃗, Size(D, D))

    κᵀ = transpose(κ)
    κθ = κ * θ
    σσᵀ = σ * transpose(σ)
    σσᵀB = σσᵀ * B
    σσᵀC = σσᵀ * C

    dA = -tr(σσᵀC) - dot(B, κθ) + 1 / 2 * dot(B, σσᵀB) - ξ₀
    dB = κᵀ * B + 2 * C * σσᵀB - 2 * C * κθ - ξ₁
    dC = κᵀ * C + C * κ + 2 * transpose(C) * σσᵀC - ξ₂

    return vcat(SVector(dA), dB, SVector(dC))
end
