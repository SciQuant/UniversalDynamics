
struct DynamicalSystem{IIP,D,M,DN,T,F,G,A,P} <: AbstractDynamics{IIP,D,M,DN,T}
    f::F
    g::G
    attributes::A
    params::P

    function DynamicalSystem{IIP,D,M,DN,T}(f::F, g::G, attrs::A, params::P) where {IIP,D,M,DN,T,F,G,A,P}
        return new{IIP,D,M,DN,T,F,G,A,P}(f, g, attrs, params)
    end
end

function DynamicalSystem(f, g, params)
    dynamics = filter(d -> d isa AbstractDynamics, values(params))
   return DynamicalSystem(f, g, dynamics, params)
end

function DynamicalSystem(f, g, dynamics, params)

    if isempty(dynamics)
        throw(ArgumentError("provide at least one `AbstractDynamics`."))
    end

    IIP = all(isinplace.(dynamics))
    OOP = all((!isinplace).(dynamics))

    if isequal(IIP, OOP)
        error("all dynamics *must* be either In-Place or Out-Of-Place.")
    end

    D = sum(dimension.(dynamics))
    M = sum(noise_dimension.(dynamics))

    DN = !any((!diagonalnoise).(dynamics))

    t0s = [initialtime(d) for d in dynamics]
    if !all(t -> t == first(t0s), t0s)
        error("all dynamics *must* have the same initial time.")
    end
    t0 = initialtime(first(dynamics))

    T = typeof(t0)

    x0 = IIP ? vcat(state.(dynamics)...) : vcat(SVector.(state.(dynamics))...)

    ρ = cat(cor.(dynamics)..., dims = (1, 2))
    ρ = IIP ? Array{T,2}(ρ) : SMatrix{size(ρ)...,T}(ρ)

    if DN
        # DiagonalNoise
        noise_rate_prototype = nothing

    elseif isone(M)
        # ScalarNoise:
        # El noise rate prototype sera una matrix de una columna y D filas (i.e. un vector
        # de dimension D). Sin embargo, en este caso el noise_rate_prototype sale de la
        # condicion inicial u0, por lo que no hace falta enviar noise_rate_prototype.
        noise_rate_prototype = nothing

    else
        prototypes = []
        for dynamic in dynamics
            if diagonalnoise(dynamic)
                # For DiagonalNoise, gprototype is a vector. But, when mixed with other
                # dynamics, it should be casted to a diagonal matrix.
                push!(prototypes, Diagonal(gprototype(dynamic)))
            else
                push!(prototypes, gprototype(dynamic))
            end
        end
        noise_rate_prototype = cat(prototypes...; dims=(1, 2))

        if IIP
            noise_rate_prototype = noise_rate_prototype # sparse(noise_rate_prototype)
        else
            noise_rate_prototype = convert(SMatrix{D,M}, noise_rate_prototype)
        end
    end

    #! Quizas la definicion del noise aqui no es tan necesaria y puede hacerse antes de
    #! simular. El problema es que en esta instancia no conozco el solver y puedo estar
    #! mandando que exista un second noise process en vano.
    noise = diffeqnoise(t0, ρ, IIP, D, M, DN)

    attrs = DynamicsAttributes(t0, x0, ρ, noise, noise_rate_prototype)

    return DynamicalSystem{IIP,D,M,DN,T}(f, g, attrs, params)
end

DynamicalSystem(dynamics) = DynamicalSystem(nothing, nothing, dynamics, nothing) # ver si mando dynamics a params tmb

for method in (:initialtime, :state, :cor, :noise, :noise_rate_prototype)
    @eval begin
        $method(ds::DynamicalSystem) = $method(ds.attributes)
    end
end
