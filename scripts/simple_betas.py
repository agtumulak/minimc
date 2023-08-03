#!/usr/bin/env python
from pyminimc import plotting, util
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

E = 0.1  # eV
kB = 8.617333262e-5  # eV / K
awr = 0.999167339
b_cutoff = 20

if __name__ == "__main__":
    # specify betas as value in [0, 1) corresponding to b_min and b_cutoff
    betas_norm = {450: [0.3]}
    Gs = [0.5]
    # append endpoints
    betas_norm = {k: [0.0] + betas_norm[k] + [1.0] for k in betas_norm}
    Gs = [0.0] + Gs + [1.0]
    betas = [
        -E / (kB * T) + norm * (b_cutoff + E / (kB * T))
        for T in betas_norm
        for norm in betas_norm[T]
    ]
    df = pd.Series(
        betas,
        index=pd.MultiIndex.from_product(
            [[E], betas_norm.keys(), Gs], names=["E", "T", "CDF"]
        ),
    ).unstack(level=["E", "T"])

    # plot
    print(df)
    x = df.iloc[:, 0].values
    y = df.index.values
    interp_x, interp_y, interp_yp = plotting.interpolate_monotonic_cubic(x, y)
    fig, ax1 = plt.subplots()
    ax1.plot(interp_x, interp_y, color="C0")
    ax1.scatter(x, y, color="C0")
    ax2 = ax1.twinx()
    ax2.plot(interp_x, interp_yp, color="C1")
    ax2.set_ylim(bottom=0)
    plt.show()

    # print corresponding outgoing energies
    energy_df = E + df * kB * df.columns.get_level_values("T")
    print("energies")
    print(energy_df)

    # preprocess for storing
    df += df.columns.get_level_values("E") / (
        kB * df.columns.get_level_values("T")
    )
    df = df.iloc[1:-1]
    df = np.log(df)

    # save to DataFrame
    U_df, S_df, V_df = util.to_svd_dfs(df)
    print(U_df)
    print(S_df)
    print(V_df)
    prefix = "/Users/atumulak/Developer/minimc/data/tnsl/debug/log_offset_beta_"
    U_df.to_hdf(prefix + "CDF.hdf5", "pandas")
    S_df.to_hdf(prefix + "S.hdf5", "pandas")
    V_df.to_hdf(prefix + "E_T.hdf5", "pandas")
