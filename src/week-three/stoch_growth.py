import dolo

from dolo.compiler.model import Model
from typing import Dict, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

model = dolo.yaml_import("src/week-three/models/sgm.yaml")

def euler_equation_error(model: Model, k_space):

    calib = model.calibration


    m = calib['exogenous']
    x = calib['controls']
    p = calib['parameters']
    ss = calib["states"]
    f = model.functions['arbitrage']

    res = []

    for k in k_space:
        s = np.array(k).reshape(ss.shape)
        res = f(m, s, x, m, s, x, p)

        np.log10(np.abs(1 - res))

    return res

k_space = np.linspace(0.01, 1, 10_000)
res = euler_equation_error(model, k_space)
plt.plot(k_space, res)
plt.savefig("test.png")

