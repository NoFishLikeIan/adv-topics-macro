import dolo

from dolo.compiler.model import Model
from typing import Dict, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def euler_equation_error(model_path: str):

    model = dolo.yaml_import(model_path)
    dr = dolo.time_iteration(model, verbose = False)

    m = model.calibration['exogenous']
    x = model.calibration['controls']
    p = model.calibration['parameters']
    ss = model.calibration["states"]
    f = model.functions['arbitrage']

    k_space = np.linspace(0.01, ss*10, num = 5_000)

    bounds = (np.min(k_space), np.max(k_space))
    N = len(k_space)

    tab = dolo.algos.simulations.tabulate(
        model, dr, "k",
        bounds=bounds,
        n_steps=N)

    c = tab["c"].tolist()

    res = []

    for i, k in enumerate(k_space):
        s = np.array(k).reshape(ss.shape)
        euler = f(m, s, x, m, s, x, p)

        rel_error = 1 - (euler + 1) / c[i]

        eee = np.log10(np.abs(rel_error))
        res.append(eee)


    return res, k_space

res, k_space = euler_equation_error(model)
plt.plot(k_space, res)
plt.savefig("test.png")

