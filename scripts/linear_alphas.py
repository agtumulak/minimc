#!/usr/bin/env python
import pandas as pd
import numpy as np
import pyminimc as pmc
from pyminimc import util
from matplotlib import pyplot as plt

E = 0.56  # eV
kB = 8.617333262e-5  # eV / K
awr = 0.999167339

beta = -0.1
T = 450

def E_p(E, beta, T):
    return E + beta * kB * T


def a_max(E, beta, T):
    return (np.sqrt(E) + np.sqrt(E + beta * kB * T)) ** 2 / (awr * kB * T)


def a_min(E, beta, T):
    return (np.sqrt(E) - np.sqrt(E + beta * kB * T)) ** 2 / (awr * kB * T)


def alpha(mu, E, beta, T):
    E_prime = E_p(E, beta, T)
    return (E + E_prime - 2 * mu * np.sqrt(E * E_prime)) / (awr * kB * T)


alpha = (E +E_p(E, beta, T))/ (awr * kB * T)
F = 1/4
alphas = pd.DataFrame(
    data=[alpha],
    columns=pd.MultiIndex.from_product(
        [[np.abs(beta)], [T]], names=("beta", "T")
    ),
    index=pd.Index([F], name="CDF"),
)

# visualize
xs = np.array([0, alpha, a_max(E, beta, T)])
ys = np.array([0, F, 1.0])
delta_xs = np.concatenate(([0], np.diff(xs)))
delta_ys = np.concatenate(([0], np.diff(ys)))
cs = 2 * delta_ys / delta_xs
cs[0] = 0
signed_cs = np.array([(-1) ** m * cs[m] for m in range(cs.size)])
cumsum_signed_cs = np.cumsum(signed_cs)
delta_xs_pow4 = delta_xs**4
numerator_factors = 2 * cumsum_signed_cs - signed_cs
numerator_terms = numerator_factors * delta_xs_pow4
f0 = -0.5 * numerator_terms.sum() / delta_xs_pow4.sum()
fs = np.array(
    [
        (-1) ** m * (f0 + cumsum_signed_cs[m])
        for m in range(len(cumsum_signed_cs))
    ]
)
plt.scatter(xs, ys)
x_points = np.linspace(xs[0], xs[-1], 100000)[:-1]
x_hi_i = np.searchsorted(xs, x_points, side="right")
x_hi = xs[x_hi_i]
x_lo = xs[x_hi_i - 1]
dydx_hi = fs[x_hi_i]
dydx_lo = fs[x_hi_i - 1]
r = (x_points - x_lo) / (x_hi - x_lo)
y_lo = ys[x_hi_i - 1]
interped = y_lo + (x_hi - x_lo) * (
    dydx_lo * r + 0.5 * (dydx_hi - dydx_lo) * r * r
)
plt.plot(x_points, interped)
plt.show()

print(alphas)

U_df, S_df, V_df = util.to_svd_dfs(alphas)
U_df.to_hdf("data/tnsl/quadratic/CDF.hdf5", "pandas")
S_df.to_hdf("data/tnsl/quadratic/S.hdf5", "pandas")
V_df.to_hdf("data/tnsl/quadratic/beta_T.hdf5", "pandas")
