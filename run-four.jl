using Plots, Roots

# Environment variables
const scriptpath = rsplit(@__FILE__, "/", limit=2)[1]
const do_plot = false
verbose = get!(ENV, "VERBOSE", "false") == "true"
plotpath = joinpath(scriptpath, "src/week-four/solutions/plots/")

include("src/week-four/partial.jl")
include("src/week-four/general.jl")

# --- Running partial equilibrium

r = .05
w = 1.
print("--> Solving partial equilibrium (r = $r, w = $w)\n")

n_steps = 25

upperbound = Float64(n_steps)

model = Aiyagari(
        0.9, 0.1, 7, # AR process parameters ρ, σ, N
        0, 0.95, 0.33, 0.1, 2 # Model parameters a, β, α, δ, σ_u
    )
    

append_filename = r != .05 ? "_$(Int(r * 100))" : ""

end_grid = true

verbose && print("Finding policy function with $(end_grid ? "endogenous grid method" : "policy function iteration")...\n")

a′, a_grid = policysolve(
    model, r, w; 
    n_steps=n_steps, upperbound=upperbound, 
    verbose=verbose, end_grid=end_grid
)


# Plot the policy function
plot(title="Policy", xaxis="a", legend=:left)

for y in model.y.transformation.(model.y.S)
    plot!(
            a_grid,
            (a -> a′(a, y)).(a_grid),
            label="a′(a | y = $(@sprintf("%.2f", y)))",
            xaxis="a"
        )
end

savefig("src/week-four/solutions/plots/partial/policy$(end_grid ? "_end" : "_pfi")$append_filename.png")


verbose && print("Finding steady state distribution...\n")

λ_eig =  distribution_pdf(a′, a_grid, model; mc=false, verbose=verbose)

xs = range(model.a_, 10., length=100)
eig_ys = λ_eig.(xs)

plot(xs, eig_ys, color=:black, label="Eig", alpha=0.5)


savefig("src/week-four/solutions/plots/partial/stat_dist$append_filename.png")


# General equilibrium

print("--> Solving general equilibrium\n")

A, K, w = solvegeneral(model, verbose=verbose, n_steps=75)

@time r = find_zero(
    x -> A(x) - K(x),
    (1 / model.β - 1  - model.δ) / 2 
)

print("Equilibrium interest rate, $r, and wage, $(w(r))\n\n")

xs = collect(range(-0.06, 0.05, length=50))
as = A.(xs)
ks = K.(xs)

plot(xs, as, label="A(r)", title="Capital market", xaxis="r")
plot!(xs, ks, label="K(r)")
savefig("src/week-four/solutions/plots/partial/k_market.png")
