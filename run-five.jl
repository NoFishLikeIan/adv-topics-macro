include("src/commons/interpolation.jl")
include("src/commons/matrix.jl")

include("src/week-five/markov/state.jl")
include("src/week-five/denhaan.jl")
include("src/week-five/policy/endgrid.jl")
include("src/week-five/krusell-smith/simulate.jl")
include("src/week-five/krusell-smith/main.jl")


using Plots

# default(dpi=600)

do_policy = true
plot_path = "src/week-five/solutions/plots/"
model = DenHaanModel()

# --- Test policy function 
if do_policy
    N_a = 100
    N_m = 10

    print("Running with default and $N_a 'a's and $N_m 'm's\n")
    
    ho_Ψ(b0, b1) = (z, K) -> exp(b0 + b1 * log(K))
    Ψ = ho_Ψ(1., 1.)
    
    g = endgrid_method(
        Ψ, model, (N_a, N_m);
        grid_bounds=[.01, 10.],
        ρ=0.5, verbose=true)

    fix_m = 1.
    as = range(0.01, 10., length=100)

    plot(title="Policy function", legend=:left, xlabel="a", ylabel="a′(a)")

    for (z_f, ϵ_f) in model.ζ.S
        ys_pol = g.(as, fix_m, z_f, collect(Float64, ϵ_f))
        plot!(as, ys_pol,
            label="a′(a | $z_f, $ϵ_f)")
    end

    savefig("$plot_path/policy.png")
end

