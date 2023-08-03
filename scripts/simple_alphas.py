#!/usr/bin/env python
from pyminimc import util
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

E = 0.1  # eV
beta = -0.1
T = 450  # K
kB = 8.617333262e-5  # eV / K
awr = 0.999167339
Ep = E + beta * kB * T  # eV


def alpha_at(mu, E, beta, T):
    Ep = E + beta * kB * T
    return (E + Ep - 2 * mu * np.sqrt(E * Ep)) / (awr * kB * T)


def mu_at(alpha, E, beta, T):
    Ep = E + beta * kB * T
    return (E + Ep - alpha * awr * kB * T) / (2 * np.sqrt(E * Ep))


if __name__ == "__main__":
    # specify alpha as mu in (-1, 1), must be in descending order
    mus = {450: [0.5]}
    Gs = [0.5]
    # append endpoints
    alphas = [alpha_at(mu, E, beta, T) for T in mus for mu in mus[T]]
    df = pd.Series(
        alphas,
        index=pd.MultiIndex.from_product(
            [[beta], mus.keys(), Gs], names=["beta", "T", "CDF"]
        ),
    ).unstack(level=["beta", "T"])

    # plot
    plottable = df.copy()
    plottable.loc[0.0] = 0.0
    plottable.loc[1.0] = 20.0  # alpha_cutoff
    plottable = plottable.sort_index()
    plottable.loc[2 * plottable.index[-1] - plottable.index[-2]] = (
        2 * plottable.iloc[-1] - plottable.iloc[-2]
    )
    plottable.loc[2 * plottable.index[0] - plottable.index[1]] = (
        2 * plottable.iloc[0] - plottable.iloc[1]
    )
    plottable = plottable.sort_index()
    # compute derivatives
    secants = (
        plottable.diff().iloc[1:].divide(np.diff(plottable.index), axis="rows")
    ) ** (-1)
    m = 3  # must satisfy 2 <= m <= 3
    s = secants.iloc[:-1]
    t = secants.shift(-1).iloc[:-1]
    u = s.where(s.abs() < t.abs(), t)
    v = s.where(s.abs() > t.abs(), t)
    dydx = np.sign(s) * m * u / (1 + (m - 1) * u / v)
    # interpoalte between coarse values
    coarse_x = plottable.iloc[:, 0].values[1:-1]
    coarse_y = plottable.index.values[1:-1]
    coarse_dydx = dydx.iloc[:, 0].values
    cubic_spline = sp.interpolate.CubicHermiteSpline(
        coarse_x, coarse_y, coarse_dydx, extrapolate=False
    )
    a_min = alpha_at(+1, E, beta, T)
    a_max = alpha_at(-1, E, beta, T)
    Hhat_min = cubic_spline(a_min)
    Hhat_max = cubic_spline(a_max)
    # plt.axvline(a_min, color="k", linestyle=":", linewidth=0.5)
    # plt.axvline(a_max, color="k", linestyle=":", linewidth=0.5)
    # plt.axhline(Hhat_min, color="k", linestyle=":", linewidth=0.5)
    # plt.axhline(Hhat_max, color="k", linestyle=":", linewidth=0.5)
    fine_alpha = np.linspace(a_min, a_max, 10000)
    fine_y = cubic_spline(fine_alpha)
    fine_dydx = cubic_spline.derivative()(fine_alpha)
    fine_palpha = cubic_spline.derivative()(fine_alpha) / (Hhat_max - Hhat_min)
    # fine_mu = mu_at(fine_alpha, E, beta, T)
    # jacobian = 2 * np.sqrt(E * Ep) / (awr * kB* T)

    # plot
    fig, [ax1, ax2] = plt.subplots(ncols=2)
    ax1.plot(fine_alpha, fine_y, color="C0")
    ax1.scatter(coarse_x, coarse_y, color="C0")
    ax2.plot(fine_alpha, fine_palpha, color="C1")
    ax2.set_ylim(bottom=0)
    plt.show()

    # print corresponding mus
    mu_df = (E + Ep - df * awr * kB * df.columns.get_level_values("T")) / (
        2 * np.sqrt(E * Ep)
    )
    print("mus")
    print(mu_df)

    # preprocess for storing
    df = np.log(df)

    # save to DataFrame
    U_df, S_df, V_df = util.to_svd_dfs(df)
    print(U_df)
    print(S_df)
    print(V_df)
    prefix = "/Users/atumulak/Developer/minimc/data/tnsl/debug/log_alpha_"
    U_df.to_hdf(prefix + "CDF.hdf5", "pandas")
    S_df.to_hdf(prefix + "S.hdf5", "pandas")
    V_df.to_hdf(prefix + "beta_T.hdf5", "pandas")
