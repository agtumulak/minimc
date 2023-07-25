#!/usr/bin/env python
import pandas as pd
import numpy as np
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


def alpha_at(mu, E, beta, T):
    E_prime = E_p(E, beta, T)
    return (E + E_prime - 2 * mu * np.sqrt(E * E_prime)) / (awr * kB * T)

def mu_at(alpha, E, beta, T):
    E_prime = E_p(E, beta, T)
    return (E + E_prime - alpha * awr * kB * T) / (2 * np.sqrt(E * E_prime))

mus = [0.0, -0.25]
alphas = [alpha_at(mu, E, beta, T) for mu in mus]
Fs = [0.5, 0.75]
alpha_df = pd.DataFrame(
    data=alphas,
    columns=pd.MultiIndex.from_product(
        [[np.abs(beta)], [T]], names=("beta", "T")
    ),
    index=pd.Index(Fs, name="CDF"),
)

# compute derivatives
xs = np.concatenate([[0.], alphas, [a_max(E, beta, T)]])
ys = np.concatenate([[0.], Fs, [1.]])
delta_xs = np.concatenate(([0], np.diff(xs)))
delta_ys = np.concatenate(([0], np.diff(ys)))
cs = 2 * delta_ys / delta_xs
cs[0] = 0
cs[1] = 0
signed_cs = np.array([(-1) ** m * cs[m] for m in range(cs.size)])
cumsum_signed_cs = np.cumsum(signed_cs)
delta_xs_pow4 = delta_xs**4
numerator_factors = 2 * cumsum_signed_cs - signed_cs
numerator_terms = numerator_factors * delta_xs_pow4
cs[1] = 0.5 * numerator_terms[2:].sum() / delta_xs_pow4[2:].sum()
fs = np.array(
    [
        (-1) ** m * (-cs[1] + cumsum_signed_cs[m])
        for m in range(len(cumsum_signed_cs))
    ]
)
fs[0] = 0
xs[0] = xs[1] - 2 * (ys[1] - ys[0]) / cs[1]
print(f"fs: {fs}")
print(f"xs: {xs}")
print(f"mus: {[mu_at(a, E, beta, T) for a in xs]}")
print(f"a_min: {a_min(E, beta, T)}")
x_points = np.linspace(xs[0], xs[-1], 1000)[:-1]
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
interped_deriv = (dydx_lo + (dydx_hi - dydx_lo) * r)

# plot derivatives
fig, ax1 = plt.subplots()
ax1.plot(x_points, interped)
ax1.scatter(xs, ys)
ax2 = ax1.twinx()
ax2.plot(x_points, interped_deriv, color='C1')
ax2.scatter(xs, fs, color='C1')
ax2.set_ylim(bottom=0)
plt.show()

# print dataframe with appended values
alpha_df_view = alpha_df.copy()
alpha_df_view.loc[0.,:] = 0.0
alpha_df_view.loc[1.,:] = a_max(E, beta, T)
print(alpha_df_view.sort_index())

U_df, S_df, V_df = util.to_svd_dfs(alpha_df)
print(U_df)
print(S_df)
print(V_df)
U_df.to_hdf("data/tnsl/quadratic/CDF.hdf5", "pandas")
S_df.to_hdf("data/tnsl/quadratic/S.hdf5", "pandas")
V_df.to_hdf("data/tnsl/quadratic/beta_T.hdf5", "pandas")
