
# From DynamicalSystemsBase.jl

printlimited(io, x::Number; Δx = 0, Δy = 0) = print(io, x)

function printlimited(io, x; Δx = 0, Δy = 0)
    sz = displaysize(io)
    io2 = IOBuffer(); ctx = IOContext(io2, :limit => true, :compact => true,
    :displaysize => (sz[1]-Δy, sz[2]-Δx))
    Base.print_array(ctx, x)
    s = String(take!(io2))
    s = replace(s[2:end], "  " => ", ")
    Base.print(io, "["*s*"]")
end

function Base.show(io::IO, ad::AbstractDynamics)
    ps = 18
    text = summary(ad)
    u0 = state(ad)'

    ctx = IOContext(io, :limit => true, :compact => true, :displaysize => (10,50))

    println(io, text)
    prefix = rpad(" state: ", ps)
    print(io, prefix); printlimited(io, u0, Δx = length(prefix)); print(io, "\n")
    println(io,  rpad(" in-place? ", ps),        isinplace(ad))
    println(io,  rpad(" Dimension: ", ps),       dimension(ad))
    println(io,  rpad(" Noise dimension: ", ps), noise_dimension(ad))
    print(io,    rpad(" diagonal noise? ", ps),  diagonalnoise(ad))
end

noise_summary(ad::AbstractDynamics) = diagonalnoise(ad) ? "DiagonalNoise" : "NonDiagonalNoise"

function Base.summary(sd::SystemDynamics)
    return "$(dimension(sd))-dimensional system dynamics with $(noise_dimension(sd))-dimensional $(noise_summary(sd))"
end

short_rate_summary(sr::OneFactorAffineModelDynamics) = "One-Factor Affine"
short_rate_summary(sr::MultiFactorAffineModelDynamics) = "Multi-Factor Affine"
short_rate_summary(sr::OneFactorQuadraticModelDynamics) = "One-Factor Quadratic"
short_rate_summary(sr::MultiFactorQuadraticModelDynamics) = "Multi-Factor Quadratic"

function Base.summary(srmd::ShortRateModelDynamics)
    return "$(dimension(srmd))-dimensional $(short_rate_summary(srmd)) Short Rate model dynamics"
end

function Base.summary(ds::DynamicalSystem)
    return "$(dimension(ds))-dimensional dynamical system with $(noise_dimension(ds))-dimensional $(noise_summary(ds))"
end
