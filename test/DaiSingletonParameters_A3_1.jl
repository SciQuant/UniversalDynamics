using NLsolve

# calibrated parameters from Dai-Singleton (2000) for the 𝔸₁(3)ₘₐₓ Ar representation:
const μ = 0.366
const ν = 0.228
const κ_rυ = 0.0348
const κ = 18.
const ῡ = 0.0158
const θ̄ = 0.0827
const σ_θυ = 0.0212
const σ_rυ = 4.2
const σ_rθ = -3.77
const σ_θr = -0.0886
const ζ2 = 0.000208
const α_r = 3.26e-14
const η2 = 0.00839
const β_θ = 7.9e-10

# useful intermediate variables
const η = sqrt(η2)
const ζ = sqrt(ζ2)
const q = (κ - ν) / κ

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

# handlers
const α₂ = x.zero[1]
const α₃ = x.zero[2]
const β₂₁ = x.zero[3]
const β₃₁ = x.zero[4]
const δ₀ = x.zero[5]
const δ₁ = x.zero[6]
const θ₁ = x.zero[7]
const κ₁₁ = x.zero[8]
const κ₂₂ = x.zero[9]
const κ₃₃ = x.zero[10]
const σ₂₁ = x.zero[11]
const σ₂₃ = x.zero[12]
const σ₃₁ = x.zero[13]
const σ₃₂ = x.zero[14]

# initial conditions (not defined in Dai Singleton paper, but needed for simulations)
# given dummy initial conditions for the AY representation:
const Y₁₀ = 1.
const Y₂₀ = 1.
const Y₃₀ = 1.
const Y₀ = [Y₁₀, Y₂₀, Y₃₀]

const L = [β₃₁ * (1 + σ₂₃)^2 0 0; 0 q 0; δ₁ 1 1]
const ϑ = [0, δ₀ + δ₁ * θ₁, δ₀]

# compute factors initial conditions for Ar representation, i.e. the vector
# given by: [υ₀, θ₀, r₀]
const υ₀, θ₀, r₀ = L * Y₀ + ϑ
