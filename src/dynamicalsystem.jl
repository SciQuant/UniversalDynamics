
# TODO: hacer esto:
#! creo que una buena idea es que esto se llame DynamicsAttributes y que este en todos los
#! AbstractDynamics que tengamos. Entonces solo 1 vez tengo que definir el isinplace, state,
#! etc! Creo que asi va a quedar mucho mas prolijo!
struct DynamicalSystemAttributes{IIP,D,M,DN,T,S,R,N,ND}
    t0::T
    x0::S
    ρ::R
    noise::N
    noise_rate_prototype::ND
end

DynamicalSystemAttributes(p::NamedTuple) = DynamicalSystemAttributes(values(p))

function DynamicalSystemAttributes(p::Tuple)
    dynamics = filter(d -> d isa AbstractDynamics, p)
    if isempty(dynamics)
        throw(ArgumentError("the container does not have any <: AbstractDynamics."))
    end
   return DynamicalSystemAttributes(dynamics)
end

function DynamicalSystemAttributes(dynamics::Tuple{Vararg{AbstractDynamics}})

    IIP = all(isinplace.(dynamics))
    OOP = all((!isinplace).(dynamics))

    # si no es ninguno de los dos, es porque estan mal dadas las condiciones inciales
    if isequal(IIP, OOP)
        error("all dynamics *must* be either `IIP` or `OOP`.")
    end

    D = sum(dimension.(dynamics))
    M = sum(noise_dimension.(dynamics))

    DN = !any((!diagonalnoise).(dynamics))

    t0s = [d.t0 for d in dynamics]
    if !all(t -> t == first(t0s), t0s)
        error("all dynamics *must* have the same initial time.")
    end
    t0 = first(dynamics).t0

    x0 = IIP ? vcat(state.(dynamics)...) : vcat(SVector.(state.(dynamics))...)

    ρ = cat(cor.(dynamics)..., dims = (1, 2))
    ρ = IIP ? ρ : SMatrix{size(ρ)...}(ρ)

    # if isequal(ρ, I)
    #     noise = nothing
    # else
    #     # IDEA: podria crear y guardar tambien el noise que no tiene condicion inicial para
    #     # el extra noise process W.dZ, i.e., el ultimo argumento le pongo nothing. De esta
    #     # forma, voy a evitar el sampling de W.dZ en los algoritmos que no necesiten el
    #     # extra process, ya que puedo mandar al solver el noise process con nothing como
    #     # condicion inicial del extra. De la forma que esta ahora, es muy probable que los
    #     # W.dZ se sorteen, pero no se use, en los algoritmos que no necesitan el extra
    #     # process. Es muy importante entender que estos se sortean antes de perform_step!,
    #     # no es que se sortea lo que se necesita en cada step.
    #     noise = IIP ?
    #         CorrelatedWienerProcess!(ρ, t0, zeros(D), zeros(D)) :
    #         CorrelatedWienerProcess(ρ, t0, @SVector(zeros(D)), @SVector(zeros(D))) # nunca probe esto, analizar

    #     # TODO: scalar noise deberia tambien construir un noise y esto deberia ir por algun
    #     # lado. O por ahi magicamente ya lo estoy handleando.
    # end

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
        # TODO: esto funciona bien cuando tenemos mezclados systemas con DN y con NDN?
        # evidentemente no, probar:
        # m1 = rand(3, 4)
        # v1 = rand(3)
        # cat(m1, v1, dims = (1, 2))
        # a aquellos que vienen de una dynamica con DN, habria que poner Diagonal(v1),
        # aunque estoy viendo que cat me devuelve un sparse, cosa que no necesariamente
        # quiero, por lo que quizas deberia hacer Matrix(Diagonal(v1))
        prototypes = []
        for dynamic in dynamics
            if diagonalnoise(dynamic)
                push!(prototypes, Diagonal(gprototype(dynamic)))
            else
                push!(prototypes, gprototype(dynamic))
            end
        end
        noise_rate_prototype = cat(prototypes...; dims=(1, 2))

        if IIP
            #! la verdad que no siempre conviene que sea sparse, por lo que deberia ser un
            #! argumento que de la opcion de hacerla o no sparse
            noise_rate_prototype = sparse(noise_rate_prototype)
        else
            noise_rate_prototype = convert(SMatrix{D,M}, noise_rate_prototype)
        end
    end

    if DN
        if isequal(ρ, I)
            noise = nothing
        else
            if IIP
                noise = CorrelatedWienerProcess!(ρ, t0, zeros(M), zeros(M))
            else
                noise = CorrelatedWienerProcess(ρ, t0, @SVector(zeros(M)), @SVector(zeros(M))) #! nunca probe con SArrays, analizar
            end
        end
    else
        if isequal(D, M)
            if isequal(ρ, I)
                # NonDiagonalNoise with a square matrix and no correlations
                noise = nothing
            else
                if IIP
                    noise = CorrelatedWienerProcess!(ρ, t0, zeros(M), zeros(M))
                else
                    noise = CorrelatedWienerProcess(ρ, t0, @SVector(zeros(M)), @SVector(zeros(M))) #! nunca probe con SArrays, analizar
                end
            end
        else
            if isequal(ρ, I)
                # NonDiagonalNoise with non-square matrix and no correlations
                if IIP
                    noise = WienerProcess!(t0, zeros(M), zeros(M))
                else
                    noise = WienerProcess(t0, @SVector(zeros(M)), @SVector(zeros(M)))
                end
            else
                if IIP
                    noise = CorrelatedWienerProcess!(ρ, t0, zeros(M), zeros(M))
                else
                    noise = CorrelatedWienerProcess(ρ, t0, @SVector(zeros(M)), @SVector(zeros(M))) #! nunca probe con SArrays, analizar
                end
            end
        end
    end

    T, S, R, N, ND = typeof.((t0, x0, ρ, noise, noise_rate_prototype))

    return DynamicalSystemAttributes{IIP,D,M,DN,T,S,R,N,ND}(
        t0, x0, ρ, noise, noise_rate_prototype
    )
end

initialtime(attrs::DynamicalSystemAttributes) = attrs.t0
state(attrs::DynamicalSystemAttributes) = attrs.x0
cor(attrs::DynamicalSystemAttributes) = attrs.ρ


struct DynamicalSystem{IIP,D,M,DN,T,F,G,A,P} <: AbstractDynamics{IIP,D,M,DN,T}
    f::F
    g::G
    attributes::A
    params::P
end

function DynamicalSystem(f, g, params)
    attrs = DynamicalSystemAttributes(params)
    #! aca se puede chequear que lo que las funciones que mandamos a dynamical system reflejen
    #! lo mismo que los atributos que se estan calculando
    return DynamicalSystem(f, g, attrs, params)
end

function DynamicalSystem(f, g, attrs::DynamicalSystemAttributes{IIP,D,M,DN,T}, params) where {IIP,D,M,DN,T}
    return DynamicalSystem{IIP,D,M,DN,T}(f, g, attrs, params)
end

# should be inner constructor?
function DynamicalSystem{IIP,D,M,DN,T}(f::F, g::G, attrs::A, params::P) where {IIP,D,M,DN,T,F,G,A,P}
    return DynamicalSystem{IIP,D,M,DN,T,F,G,A,P}(f, g, attrs, params)
end