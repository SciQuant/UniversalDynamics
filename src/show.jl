
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

function Base.show(io::IO, ds::AbstractDynamics)
    ps = 18
    text = summary(ds)
    u0 = state(ds)'

    ctx = IOContext(io, :limit => true, :compact => true, :displaysize => (10,50))

    println(io, text)
    prefix = rpad(" state: ", ps)
    print(io, prefix); printlimited(io, u0, Δx = length(prefix)); print(io, "\n")
    println(io,  rpad(" in-place? ", ps),        isinplace(ds))
    println(io,  rpad(" Dimension: ", ps),       dimension(ds))
    println(io,  rpad(" Noise dimension: ", ps), noise_dimension(ds))
    print(io,    rpad(" diagonal noise? ", ps),  diagonalnoise(ds))
end

function Base.summary(sd::SystemDynamics)
    noise = diagonalnoise(sd) ? "DiagonalNoise" : "NonDiagonalNoise"
    return "$(dimension(sd))-dimensional system dynamics with $(noise_dimension(sd))-dimensional $(noise)"
end

function Base.summary(ds::DynamicalSystem)
    noise = diagonalnoise(ds) ? "DiagonalNoise" : "NonDiagonalNoise"
    return "$(dimension(ds))-dimensional dynamical system with $(noise_dimension(ds))-dimensional $(noise)"
end
