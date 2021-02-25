
import StochasticDiffEq: SDEProblem, solve

function SDEProblem(ds::DynamicalSystem{IIP}, tspan; u0=get_state(ds), kwargs...) where {IIP}
    return SDEProblem{IIP}(
        SDEFunction(ds.f, ds.g),
        ds.g,
        u0,
        tspan,
        get_parameters(ds);
        # noise=get_noise(ds),
        noise_rate_prototype=get_noise_rate_prototype(ds),
        kwargs...
    )
end

function solve(ds::DynamicalSystem, T; u0=get_state(ds), alg=SRIW1(), diffeq...)
    noise = diffeqnoise(ds, alg)
    prob = SDEProblem(ds, (get_t0(ds), T); u0=u0, noise=noise)
    sol = solve(prob, alg; diffeq...)
    return sol
end

function montecarlo(
    ds::DynamicalSystem, T, trials::Integer=1;
    u0=get_state(ds), alg=SRIW1(), ensemblealg=EnsembleThreads(), seed::Integer=314159, diffeq...
)
    tspan = (get_t0(ds), T)
    noise = diffeqnoise(ds, alg)
    prob = SDEProblem(ds, tspan; u0=u0, noise=noise)

    rng = Xorshifts.Xoroshiro128Plus(seed)
    seeds = rand(rng, UInt64, trials)

    prob_func = (prob, i, repeat) -> begin
        return StochasticDiffEq.remake(prob, seed = seeds[i])
        # remake(prob, p = copy.(prob.p), seed = seeds[i]) # thread safe (?)
    end

    sde = EnsembleProblem(prob, prob_func=prob_func, safetycopy=true)
    sol = solve(sde, alg, ensemblealg; trajectories=trials, diffeq...)

    return sol
end