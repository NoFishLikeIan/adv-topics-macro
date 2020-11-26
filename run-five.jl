include("src/commons/interpolation.jl")
include("src/commons/matrix.jl")

include("src/week-five/markov/state.jl")
include("src/week-five/denhaan.jl")
include("src/week-five/policy/pfi.jl")
include("src/week-five/krusell-smith/main.jl")


using Plots

default(size=(12, 8), dpi=600)

do_policy = true
plot_path = "src/week-five/solutions/plots/"
model = DenHaanModel()

# --- Test policy function 
if do_policy
    N = 10
    
    ho_Ψ(b0, b1) = (z, K) -> exp(b0 + b1 * log(K))
    Ψ = ho_Ψ(0., 1.)
    
    policy = pfi(
        makeproduction(model),
        makeutility(model),
        Ψ, model, (N, N);
        ρ=0.7)

    fix_m = 1.
    as = range(0.01, 10., length=100)
    plot(
        title="Policy function", legend=:left,
        xlabel="a", ylabel="a′(a)"
    )

    for (z_f, ϵ_f) in mode.ζ.S
        ys_pol = policy.(as, fix_m, z_f, ϵ_f)
        plot!(
            as, ys_pol,
            label="a′(a | $z_f, $ϵ_f)"
        )
    end

    savefig("$plot_path/policy.png")
end
