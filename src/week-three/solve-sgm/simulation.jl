function simulate_shock(c_p, k_p, k1, shock)
    y = exp.(shock)
    T = length(y)
    
    c1 = c_p(k1, y[1])

    tab = zeros((T, 2))
    tab[1, :] = [c1, k1]

    for t in 2:T
        k = tab[t - 1, 2]
        tab[t, :] = [c_p(k, y[t]), k_p(k, y[t])]
    end

    return tab
end