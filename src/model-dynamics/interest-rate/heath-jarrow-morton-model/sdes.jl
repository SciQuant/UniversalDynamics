
@inline function drift!(du, u, p::HJMP{true,D,D,true,T,Terminal}, t) where {D,T}
    @unpack Tenors, τ, ρ, cache = p
    @unpack σ = cache
    # p.σ(σ, t, ?)
    
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
                du[i] = - σ(t, Tenors[i])[i-1] * quadgk(u -> σ(t, u)[i-1], t, Tenors[i])[1]
            end
        end
    end

    return convert(SVector, du)
end

@inline function drift!(du, u, p::HJMP{true,D,D,true,T,Spot}, t) where {D,T}
    @unpack Tenors, τ, ρ, cache = p
    @unpack σ = cache
    # p.σ(σ, t, ?)
    
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
    # p.σ(σ, t, ?)
   
    return nothing
end

@inline function diffusion(u, p::HJMP{false,D,D,true}, t) where {D}
    @unpack Tenors, τ, ρ = p
    σ = p.σ

    du = @MVector zeros(eltype(u), D)

    @inbounds begin
        for i in 2:D
            du[i] = t > Tenors[i] ? zero(eltype(du)) : σ(t, Tenors[i])[i-1]
        end
    end

    return convert(SVector, du)
end
