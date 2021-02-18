
@inline function drift!(du, u, p::HJMP{true,D,D,true,T,Terminal}, t) where {D,T}
    @unpack Tenors, τ, ρ, cache = p
    @unpack σ, σI = cache

    @inbounds begin
        # 'i' ranges from 2 to D because f₁ has μ = 0
        for i in 2:D
            du[i] = zero(eltype(u))
            if t ≤ Tenors[i]

                p.σ(σ, t, Tenors[i])
                du[i] += - σ[i - 1] 

                quadgk!((z, w) -> p.σ(z, t, w), σI, t, Tenors[i])
                du[i] *= σI[i - 1]
            end
        end
    end
    
    return nothing
end

@inline function drift(u, p::HJMP{false,D,D,true,T,Terminal}, t) where {D,T}
    @unpack Tenors, τ, ρ = p
    σ = p.σ

    du = @MVector zeros(eltype(u), D)

    @inbounds begin
        # 'i' ranges from 2 to D because f₁ has μ = 0
        for i in 2:D
            if t ≤ Tenors[i]
                du[i] = - σ(t, Tenors[i])[i - 1] * quadgk(u -> σ(t, u)[i - 1], t, Tenors[i])[1]
            end
        end
    end

    return convert(SVector, du)
end

@inline function drift!(du, u, p::HJMP{true,D,D,true,T,Spot}, t) where {D,T}
    @unpack Tenors, τ, ρ, cache = p
    @unpack σ, σI = cache

    @inbounds begin
        # 'i' ranges from 2 to D because f₁ has μ = 0
        for i in 2:D
            du[i] = zero(eltype(u))
            if t ≤ Tenors[i]

                p.σ(σ, t, Tenors[i])
                du[i] += σ[i - 1] 

                quadgk!((z, w) -> p.σ(z, t, w), σI, t, Tenors[i])
                du[i] *= σI[i - 1]
            end
        end
    end
    
    return nothing
end

@inline function drift(u, p::HJMP{false,D,D,true,T,Spot}, t) where {D,T}
    @unpack Tenors, τ, ρ = p
    σ = p.σ

    du = @MVector zeros(eltype(u), D)

    @inbounds begin
        # 'i' ranges from 2 to D because f₁ has μ = 0
        for i in 2:D
            if t ≤ Tenors[i]
                du[i] = σ(t, Tenors[i])[i-1] * quadgk(u -> σ(t, u)[i-1], t, Tenors[i])[1]
            end
        end
    end

    return convert(SVector, du)
end

@inline function diffusion!(du, u, p::HJMP{true,D,D,true}, t) where {D}
    @unpack Tenors, cache = p
    @unpack σ = cache
    
    @inbounds begin
        for i in 2:D
            p.σ(σ, t, Tenors[i])
            du[i] = t > Tenors[i] ? zero(eltype(u)) : σ[i - 1]
        end
    end
    return nothing
end

@inline function diffusion(u, p::HJMP{false,D,D,true}, t) where {D}
    @unpack Tenors, τ, ρ = p
    σ = p.σ

    du = @MVector zeros(eltype(u), D)

    @inbounds begin
        for i in 2:D
            du[i] = t > Tenors[i] ? zero(eltype(u)) : σ(t, Tenors[i])[i - 1]
        end
    end

    return convert(SVector, du)
end
