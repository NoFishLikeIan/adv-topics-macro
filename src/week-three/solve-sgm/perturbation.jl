using Dolo, Dolang, StaticArrays

function pert_eee(model; n_steps=1_000, bounds=[.01, 1.])
    res = time_iteration(model; verbose=false)

    m = model.calibration[:exogenous]
    x = model.calibration[:controls]
    p = model.calibration[:parameters]
    ss = model.calibration[:states]

    f = eval(Dolang.gen_generated_gufun(model.factories[:arbitrage]))

    tab = tabulate(model, res.dr, :k, bounds, ss, m; n_steps=n_steps)
    c = tab[:c]
    k_space = tab[:k]

    eee = zeros(size(k_space))

    for (i, k) in enumerate(k_space)

        euler = f(m, [k], x, m, [k], x, p)

        rel_error = 1 - (euler[1] + 1) / c[i]

        eee[i] = log10(abs(rel_error))

    end

    return eee, tab

end