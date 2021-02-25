
for method in (:drift, :diffusion)
    method! = Symbol(method, !)
    @eval begin
        $method(x, srmd::ShortRateModelDynamics, t) = $method(x, get_parameters(srmd), t)
        $method!(dx, x, srmd::ShortRateModelDynamics, t) = $method!(dx, x, get_parameters(srmd), t)
    end
end

############################# One Factor Affine + Out of Place ##############################

function drift(r, p::AffineParameters{OneFactor,false,D,DN}, t) where {D,DN}
    @unpack κ, θ = p
    return κ(t) * (θ(t) - r[1])
end

function diffusion(r, p::AffineParameters{OneFactor,false,D,DN}, t) where {D,DN}
    @unpack Σ, α, β = p
    return Σ(t) * sqrt(α(t) + β(t) * r[1])
end

############################### One Factor Affine + In Place ################################
#! en las funciones OneFactor con IIP = true, uso o no uso cache? Pido funciones f(u, t) o
#! f(t)? Por ahora me quedo con la segunda opcion, es mas facil!

function drift!(dr, r, p::AffineParameters{OneFactor,true,D,DN}, t) where {D,DN}
    @unpack κ, θ = p
    dr[1] = κ(t) * (θ(t) - r[1])
    return nothing
end

function diffusion!(dr, r, p::AffineParameters{OneFactor,true,D,DN}, t) where {D,DN}
    @unpack Σ, α, β = p
    dr[1] = Σ(t) * sqrt(α(t) + β(t) * r[1])
    return nothing
end

############################ Multi Factor Affine + Out of Place #############################

function drift(x, p::AffineParameters{MultiFactor,false,D,DN}, t::Real) where {D,DN}
    @unpack κ, θ = p
    return κ(t) * (θ(t) - x)
end

# diagonal noise
function diffusion(x, p::AffineParameters{MultiFactor,false,D,true}, t) where {D}
    @unpack Σ, α, β = p

    S = α(t) + β(t) * x
    σ = Σ(t) .* sqrt.(S) # Σ is a vector

    return σ
end

# non-diagonal noise
function diffusion(x, p::AffineParameters{MultiFactor,false,D,false}, t) where {D}
    @unpack Σ, α, β = p

    S = Diagonal(α(t) + β(t) * x)
    σ = Σ(t) * sqrt(S) # Σ is a matrix

    return σ
end

############################## Multi Factor Affine + In Place ###############################

function drift!(dx, x, p::AffineParameters{MultiFactor,true,D,DN}, t) where {D,DN}
    @unpack cache = p
    @unpack κ, θ, v1 = cache

    p.κ(κ, t)
    p.θ(θ, t)

    v1 .= θ .- x
    mul!(dx, κ, v1)

    return nothing
end

# diagonal noise
function diffusion!(dx, x, p::AffineParameters{MultiFactor,true,D,true}, t) where {D}
    @unpack cache = p
    @unpack Σ, α, β, v1 = cache

    p.Σ(Σ, t) # Σ is a vector
    p.α(α, t)
    p.β(β, t)

    # mul!(v1, β, x) .+= α
    copyto!(v1, α)
    mul!(v1, β, x, true, true)
    map!(sqrt, v1, v1)
    dx .= Σ .* v1

    return nothing
end

# non-diagonal noise
function diffusion!(dx, x, p::AffineParameters{MultiFactor,true,D,false}, t) where {D}
    @unpack cache = p
    @unpack Σ, α, β, v1 = cache

    p.Σ(Σ, t) # Σ is a matrix
    p.α(α, t)
    p.β(β, t)

    # mul!(v1, β, x) .+= α
    copyto!(v1, α)
    mul!(v1, β, x, true, true)
    map!(sqrt, v1, v1)
    mul!(dx, Σ, Diagonal(v1)) # permutedims is called and it allocates

    return nothing
end

########################### One Factor Quadratic + Out of Place ############################

function drift(x, p::QuadraticParameters{OneFactor,false,D,DN}, t) where {D,DN}
    @unpack κ, θ = p
    return κ(t) * (θ(t) - x)
end

function diffusion(x, p::QuadraticParameters{OneFactor,false,D,DN}, t) where {D,DN}
    return p.σ(t)
end

############################# One Factor Quadratic + In Place ##############################
#! en las funciones OneFactor con IIP = true, uso o no uso cache? Pido funciones f(u, t) o
#! f(t)? Por ahora me quedo con la segunda opcion, es mas facil!

function drift!(dx, x, p::QuadraticParameters{OneFactor,true,D,DN}, t) where {D,DN}
    @unpack κ, θ = p
    dx[1] = κ(t) * (θ(t) - x)
    return nothing
end

function diffusion!(dx, x, p::QuadraticParameters{OneFactor,true,D,DN}, t) where {D,DN}
    @unpack σ = p
    dx[1] = σ(t)
    return nothing
end

########################### Multi Factor Quadratic + Out of Place ##########################

function drift(x, p::QuadraticParameters{MultiFactor,false,D,DN}, t) where {D,DN}
    @unpack κ, θ = p
    return κ(t) * (θ(t) - x)
end

# diagonal noise
function diffusion(x, p::QuadraticParameters{MultiFactor,false,D,true}, t) where {D}
    @unpack σ = p
    return σ(t) # σ is a vector
end

# non-diagonal noise
function diffusion(x, p::QuadraticParameters{MultiFactor,false,D,false}, t) where {D}
    @unpack σ = p
    return σ(t) # σ is a matrix
end

############################# Multi Factor Quadratic + In Place ############################

function drift!(dx, x, p::QuadraticParameters{MultiFactor,true,D,DN}, t) where {D,DN}
    @unpack cache = p
    @unpack κ, θ, v1 = cache

    p.κ(κ, t)
    p.θ(θ, t)

    v1 .= θ .- x
    mul!(dx, κ, v1)

    return nothing
end

# diagonal noise
function diffusion!(dx, x, p::QuadraticParameters{MultiFactor,true,D,true}, t) where {D}
    @unpack σ = p
    copyto!(dx, σ(t)) # σ is a vector
    return nothing
end

# non-diagonal noise
function diffusion!(dx, x, p::QuadraticParameters{MultiFactor,true,D,false}, t) where {D}
    @unpack σ = p
    copyto!(σx, σ(t)) # σ is a matrix
    return nothing
end
