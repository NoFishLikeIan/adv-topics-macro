include("src/commons/interpolation.jl")
include("src/commons/matrix.jl")

include("src/week-five/markov/state.jl")
include("src/week-five/denhaan.jl")

include("src/week-five/policy/endgrid.jl")
include("src/week-five/policy/vfi.jl")

include("src/week-five/krusell-smith/simulate.jl")
include("src/week-five/krusell-smith/det_simulation.jl")
include("src/week-five/krusell-smith/main.jl")

Random.seed!(11148706)

using Plots, KernelDensity

default(dpi=300)

plot_path = "src/week-five/solutions/plots/"

# Test algorithm

function run_krusselsmith(append="", μ=0.15, plot=false)
    model = DenHaanModel(μ)
    g = krusellsmith(model; verbose=true, ρ=0.5)

    a_grid = range(0.01, 5., length=500)

    as, zs = economysim(g, model; drop=0, T=5_000)

    if plot
        plot(title="Policy function", legend=:left, xlabel="a", ylabel="a′(a)")

        for (z_f, ϵ_f) in model.ζ.S
            ys_pol = g.(a_grid, 1., z_f, collect(Float64, ϵ_f))
            plot!(a_grid, ys_pol, label="a′(a | $z_f, $ϵ_f)")
        end

        savefig("$plot_path/policy_krusell$append.png")

        boom = zs .== 1.01

        boom_as = vec(mean(as[boom, :], dims=1))
        recession_as = vec(mean(as[.!boom, :], dims=1))

        a_grid = range(0.01, 5., length=500)
        plot(title="Steady state distribution", xaxis="a", yaxis="probability")

        plot!(a_grid, pdf(kde(boom_as), a_grid), label="boom")

        plot!(a_grid, pdf(kde(recession_as), a_grid), label="recession")

        savefig("$plot_path/stationary_distribution$append.png")
    end

    return model, g, as, zs
end

function gif_plot(as, zs; append="")

    densities = []
    
    for (t, z) in enumerate(zs[end - 100:end])
        boom = z == 1.01
        color = boom ? :green : :red
    
        density = pdf(kde(as[t, :]), a_grid)
    
        density = density ./ sum(density)
    
        push!(densities, (color, density))
    
    end
    
    y_upper = maximum(maximum(den) for (c, den) in densities)
    
    anim = @animate for (t, (color, dens)) in enumerate(densities)
    
        plot(
            a_grid, dens, 
            label="$t / $(length(densities))",
            ylims=(0., y_upper),
            color=color, title="Assets distribution",
            xaxis="a", yaxis="probability"
        )
    
    end
    
    gif(anim, "$plot_path/dist$append.gif", fps=4)

end


model, g, as, zs = run_krusselsmith()
gif_plot(as, zs)