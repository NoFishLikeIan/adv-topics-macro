
"""
Simulate the economy for N individuals and T periods
    given a policy function
        
TODO: move to only one loop combining the conditional simulation 
"""
function economysim(
    g::Function, model::DenHaanModel;
    N::Int=10_000, T::Int=1_500, verbose=false)
    
    zs = simulation(model.Z, T)
    ϵs = conditional_simulation(model, zs, N)

    as = similar(ϵs)
    as[1, :] = rand(Uniform(), 1, N)

    for (t, z) in enumerate(zs[1:end - 1])
        verbose && print("Simulating economy $t / $(T - 1)\r")
        ϵ = ϵs[t, :]
        a = as[t, :]

        m = mean(a)

        a′ = @. g(a, m, z, ϵ)
        as[t + 1, :] = a′
    end

    verbose && print("\n")

    return as, zs
end

