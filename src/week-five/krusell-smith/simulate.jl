
"""
Simulate the economy for N individuals and T periods
    given a policy function
"""
function economysim(
    g::Function, model::DenHaanModel;
    N::Int=10_000, T::Int=1_500, verbose=false)

    a0 = rand(Uniform(10., 20.), N)

    ϵs = simulation(model.Ε, T)
    zs = conditional_simulation(model.ζ, ϵs, N)

    as = Array{Float64}(undef, T, N)
    as[1, :] = a0

    for (t, ϵ) in enumerate(collect(Float64, ϵs))
        verbose && print("Simulating economy $t / $T\r")
        z = zs[t, :]
        a = as[t, :]

        m = mean(a)

        a′ = @. g(a, m, z, ϵ)
        as[t, :] = a′
    end

    verbose && print("\n")

    return as

end

