#!/usr/bin/env python

import numpy as np
import pandas as pd
from pyminimc import compression, plotting
import matplotlib.pyplot as plt

beta_cutoff = 20
kB = 8.617333262e-11

if __name__ == "__main__":
    df = pd.DataFrame(
        pd.read_hdf(
            "/Users/atumulak/Developer/minimc/data/tnsl/endfb8/dfs/beta_full.hdf5"
        )
    )

    # convert to offsets from b_min
    b_min = -df.columns.get_level_values("E") / (
        kB * df.columns.get_level_values("T")
    )
    df = df - b_min

    # coarsen energies
    fine_E_grid = df.columns.get_level_values("E").unique()
    coarse_E_grid = fine_E_grid[
        np.linspace(0, fine_E_grid.size - 1, 50).astype(int)
    ]

    # coarsen temperatures
    fine_T_grid = df.columns.get_level_values("T").unique()
    coarse_T_grid = fine_T_grid

    # coarsen CDFs
    fine_CDF_grid = df.index[1:-1]  # don't include zero and one
    coarse_CDF_grid = fine_CDF_grid[
        np.linspace(0, fine_CDF_grid.size - 1, 20).astype(int)
    ]

    # select coarse indices
    coarse = df.loc[
        coarse_CDF_grid,
        pd.MultiIndex.from_product([coarse_E_grid, coarse_T_grid]),
    ]

    # plot some curves
    truncated_plottable = np.exp(compression.truncate(np.log(coarse), rank=37))
    truncated_plottable.loc[0.0] = 0.0
    truncated_plottable.loc[
        1.0
    ] = beta_cutoff + coarse.columns.get_level_values("E") / (
        kB * coarse.columns.get_level_values("T")
    )
    truncated_plottable = truncated_plottable.sort_index()
    nonmonotonic_cols = (truncated_plottable.diff().iloc[1:] < 0).any().sum()
    print(f"nonmonotonic columns: {nonmonotonic_cols}")
    for plot_E in coarse_E_grid[
        np.linspace(0, coarse_E_grid.size - 1, 25).astype(int)
    ]:
        for plot_T in coarse_T_grid[
            np.linspace(0, coarse_T_grid.size - 1, 3).astype(int)
        ]:
            truncated_col = truncated_plottable.loc[:, (plot_E, plot_T)]
            (
                interp_x,
                interp_y,
                interp_yp,
            ) = plotting.interpolate_monotonic_cubic(
                truncated_col.values, truncated_col.index, N=100000, m=3
            )
            true_col = df.loc[:, (plot_E, plot_T)]
            fig, ax1 = plt.subplots()
            ax1.plot(interp_x, interp_y, color="C0")
            ax1.scatter(truncated_col.values, truncated_col.index, color="C0")
            ax1.plot(true_col.values, true_col.index, color="C0", linestyle=":")
            ax1.set_xlim(left=-10.0, right=25)
            ax2 = ax1.twinx()
            ax2.plot(interp_x, interp_yp, color="C1")
            ax2.set_ylim(bottom=0)
            plt.title(f"E: {plot_E:.5e}, T: {plot_T}")
            plt.show()

    U_df, S_df, V_df = util.to_svd_dfs(np.log(coarse), order=37)
    reconstructed = np.exp(util.from_svd_dfs(U_df, S_df, V_df))
    U_df.to_hdf(
        "~/Developer/minimc/data/tnsl/partitionless/log_offset_beta_CDF.hdf5",
        "pandas",
    )
    S_df.to_hdf(
        "~/Developer/minimc/data/tnsl/partitionless/log_offset_beta_S.hdf5",
        "pandas",
    )
    V_df.to_hdf(
        "~/Developer/minimc/data/tnsl/partitionless/log_offset_beta_E_T.hdf5",
        "pandas",
    )

    # check for nonmonoticity as a function of rank
    for rank in range(np.linalg.matrix_rank(coarse) + 1):
        truncated = np.exp(compression.truncate(np.log(coarse), rank=rank))
        nonmonotonic = (truncated.diff() < 0).sum().sum()
        print(f"rank={rank}, nonmonotonic={nonmonotonic}")
