using Parameters, Interpolations

function endgrid(
    ai::Aiyagari, r::Float64, w::Float64;
    n_steps=1_000, upperbound=2)

    u, u_c = make_u(ai)
    @unpack β, a_, y = ai 
    Γ = y.P

    a_p = collect(range(a_, upperbound, length=n_steps))

end