## Dynamics

So far we have only mentioned the existence of [`UniversalDynamics.InterestRateModelDynamics`](@ref). We now focus on Term Structure Models, given by the following abstract type:

```@docs
UniversalDynamics.TermStructureModelDynamics
```

The library implements the following Term Structure Models:

```@docs
UniversalDynamics.ShortRateModelDynamics
```

See [`Short Rate Models`](@ref ShortRateModelIntroduction) section for detailed information.

```@docs
UniversalDynamics.LiborMarketModelDynamics
```

See [`Libor Market Model`](@ref LiborMarketModelIntroduction) section for detailed information.
