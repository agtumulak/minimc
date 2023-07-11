#!/usr/bin/env python
import pandas as pd
import numpy as np
import pyminimc as pmc
from pyminimc import util

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


# the alpha dataset that would be used
alphas = pd.DataFrame(
    data=[0.5 * a_max(E, beta,T)],
    columns=pd.MultiIndex.from_product(
        [[np.abs(beta)], [T]], names=("beta", "T")
    ),
    index=pd.Index([0.5], name="CDF"),
)

print(alphas)
print(" ")
print(f"amax: {a_max(E, beta, T)}")

U_df, S_df, V_df = util.to_svd_dfs(alphas)

U_df.to_hdf("data/tnsl/quadratic/CDF.hdf5", "pandas")
S_df.to_hdf("data/tnsl/quadratic/S.hdf5", "pandas")
V_df.to_hdf("data/tnsl/quadratic/beta_T.hdf5", "pandas")
