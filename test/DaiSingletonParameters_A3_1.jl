using NLsolve

function DaiSingletonParameters()

    # calibrated parameters from Dai-Singleton (2000) for the 𝔸₁(3)ₘₐₓ Ar representation:
    μ = 0.366
    ν = 0.228
    κ_rυ = 0.0348
    κ = 18.
    ῡ = 0.0158
    θ̄ = 0.0827
    σ_θυ = 0.0212
    σ_rυ = 4.2
    σ_rθ = -3.77
    σ_θr = -0.0886
    ζ2 = 0.000208
    α_r = 3.26e-14
    η2 = 0.00839
    β_θ = 7.9e-10

    # useful intermediate variables
    η = sqrt(η2)
    ζ = sqrt(ζ2)
    q = (κ - ν) / κ

    # compute the AY representation parameters by means of a non linear solver
    function f!(F, x)
        F[1] = μ - x[8]
        F[2] = ν - x[9]
        F[3] = κ - x[10]
        F[4] = η2 - x[4] * (1 + x[12])^2
        F[5] = ζ2 - q^2 * x[1]
        F[6] = κ_rυ - x[6] * (x[8] - x[10]) / x[4] / ((1 + x[12])^2)
        F[7] = σ_rυ - (x[6] + x[11] + x[13]) / x[4] / ((1 + x[12])^2)
        F[8] = β_θ - q^2 * x[3] / x[4] / ((1 + x[12])^2)
        F[9] = ῡ - x[4] * x[7] * (1 + x[12])^2
        F[10] = θ̄ - x[5] + x[6] * x[7]
        F[11] = σ_rθ - (1 + x[14]) / q
        F[12] = σ_θυ - q * x[11] / x[4] / ((1 + x[12])^2)
        F[13] = σ_θr - q * x[12] / (1 + x[12])
        F[14] = α_r - x[2] * (1 + x[12])^2
    end

    # solve
    x = nlsolve(f!, ones(14), autodiff = :forward)

    α₂ = x.zero[1]
    α₃ = x.zero[2]
    β₂₁ = x.zero[3]
    β₃₁ = x.zero[4]
    δ₀ = x.zero[5]
    δ₁ = x.zero[6]
    θ₁ = x.zero[7]
    κ₁₁ = x.zero[8]
    κ₂₂ = x.zero[9]
    κ₃₃ = x.zero[10]
    σ₂₁ = x.zero[11]
    σ₂₃ = x.zero[12]
    σ₃₁ = x.zero[13]
    σ₃₂ = x.zero[14]

    # initial conditions (not defined in Dai Singleton paper, but needed for simulations)
    # given dummy initial conditions for the AY representation:
    Y₁₀ = 1.
    Y₂₀ = 1.
    Y₃₀ = 1.
    Y₀ = [Y₁₀, Y₂₀, Y₃₀]

    L = [β₃₁ * (1 + σ₂₃)^2 0 0; 0 q 0; δ₁ 1 1]
    ϑ = [0, δ₀ + δ₁ * θ₁, δ₀]

    # compute factors initial conditions for Ar representation, i.e. the vector given by:
    # [υ₀, θ₀, r₀]
    υ₀, θ₀, r₀ = L * Y₀ + ϑ

    return (υ₀, θ₀, r₀, μ, ν, κ_rυ, κ, ῡ, θ̄, η, σ_θυ, σ_θr, σ_rυ, σ_rθ, ζ, α_r, β_θ)
end