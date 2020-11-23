using Plots

# Environment variables
const scriptpath = rsplit(@__FILE__, "/", limit=2)[1]
const do_plot = false
verbose = get!(ENV, "VERBOSE", "false") == "true"
plotpath = joinpath(scriptpath, "src/week-four/solutions/plots/")

include("src/week-four/partial.jl")

# --- Running partial equilibrium


r = .05
w = 1.

n_steps = 500
upperbound = 8.

model = Aiyagari(
        0.9, 0.1, 7, # AR process parameters ρ, σ, N
        0, 0.95, 0.33, 0.1, 2 # Model parameters a, β, α, δ, σ_u
    )
    

verbose && print("Finding policy function...\n")
a′, a_grid = policysolve(
    model, r, w; 
    n_steps=n_steps, upperbound=upperbound, 
    verbose=verbose
)

# Plot the policy function
plot(title="Policy", xaxis="a")

for y in model.y.transformation.(model.y.S)
    plot!(
            a_grid,
            (a -> a′(a, y)).(a_grid),
            label="a′(a | y = $(@sprintf("%.2f", y)))",
            xaxis="a"
        )
end

append_filename = r != .05 ? "_$(Int(r * 100))" : ""

savefig("src/week-four/solutions/plots/partial/policy$append_filename.png")


verbose && print("Finding steady state distribution...\n")

λ_mc = distribution_pdf(a′, a_grid, model; mc=true, verbose=verbose)
λ_eig =  distribution_pdf(a′, a_grid, model; mc=false, verbose=verbose)

mc_ys = λ_mc.(a_grid) / sum(λ_mc.(a_grid))

plot(
    a_grid, mc_ys, color=:red, label="MC", alpha=0.5,
    title="Stable distribution", xaxis="a", yaxis="λ(a)"
)

eig_ys = λ_eig.(a_grid) / sum(λ_eig.(a_grid))

plot!(a_grid, λ_eig, color=:black, label="Eig", alpha=0.5)


savefig("src/week-four/solutions/plots/partial/stat_dist$append_filename.png")
