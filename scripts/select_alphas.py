#!/usr/bin/env python

import numpy as np
import pandas as pd
from pyminimc import util, compression, plotting
from matplotlib import pyplot as plt

alpha_cutoff = 1636.7475317348378

if __name__ == "__main__":
    df = pd.DataFrame(
        pd.read_hdf(
            "/Users/atumulak/Developer/minimc/data/tnsl/endfb8/dfs/alpha_full.hdf5"
        )
    )

    # coarsen betas
    fine_beta_grid = df.columns.get_level_values("beta").unique()
    coarse_beta_grid = fine_beta_grid[
        np.linspace(0, fine_beta_grid.size - 1, 25).astype(int)
    ]

    # coarsen temperatures
    fine_T_grid = df.columns.get_level_values("T").unique()
    coarse_T_grid = fine_T_grid

    # coarsen CDFs
    fine_CDF_grid = df.index[1:-1]  # don't include zero and one
    coarse_CDF_grid = fine_CDF_grid[
        np.linspace(0, fine_CDF_grid.size - 1, 10).astype(int)
    ]

    # select coarse indices
    coarse = df.loc[
        coarse_CDF_grid,
        pd.MultiIndex.from_product([coarse_beta_grid, coarse_T_grid]),
    ]

    # plot some curves
    truncated_plottable = np.exp(compression.truncate(np.log(coarse), rank=6))
    truncated_plottable.loc[0.0] = 0.0
    truncated_plottable.loc[1.0] = alpha_cutoff
    truncated_plottable = truncated_plottable.sort_index()
    nonmonotonic_cols = (truncated_plottable.diff().iloc[1:] < 0).any().sum()
    print(f"nonmonotonic columns: {nonmonotonic_cols}")
    for plot_beta in coarse_beta_grid[
        np.linspace(0, coarse_beta_grid.size - 1, 25).astype(int)
    ]:
        for plot_T in coarse_T_grid[
            np.linspace(0, coarse_T_grid.size - 1, 3).astype(int)
        ]:
            truncated_col = truncated_plottable.loc[:, (plot_beta, plot_T)]
            (
                interp_x,
                interp_y,
                interp_yp,
            ) = plotting.interpolate_monotonic_cubic(
                truncated_col.values, truncated_col.index, N=100000, m=3
            )
            true_col = df.loc[:, (plot_beta, plot_T)]
            fig, ax1 = plt.subplots()
            ax1.plot(interp_x, interp_y, color="C0")
            ax1.scatter(truncated_col.values, truncated_col.index, color="C0")
            ax1.plot(true_col.values, true_col.index, color="C0", linestyle=":")
            ax1.set_xlim(left=-10.0, right=coarse.max().max())
            ax2 = ax1.twinx()
            ax2.plot(interp_x, interp_yp, color="C1")
            ax2.set_ylim(bottom=0)
            plt.title(f"beta: {plot_beta:.5e}, T: {plot_T}")
            plt.show()

    U_df, S_df, V_df = util.to_svd_dfs(np.log(coarse), order=4)
    reconstructed = np.exp(util.from_svd_dfs(U_df, S_df, V_df))
    U_df.to_hdf(
        "~/Developer/minimc/data/tnsl/endfb8/log_alpha_CDF_coarse.hdf5",
        "pandas",
    )
    S_df.to_hdf(
        "~/Developer/minimc/data/tnsl/endfb8/log_alpha_S_coarse.hdf5",
        "pandas",
    )
    V_df.to_hdf(
        "~/Developer/minimc/data/tnsl/endfb8/log_alpha_beta_T_coarse.hdf5",
        "pandas",
    )

    # check for nonmonoticity as a function of rank
    for rank in range(np.linalg.matrix_rank(coarse) + 1):
        truncated = np.exp(compression.truncate(np.log(coarse), rank=rank))
        nonmonotonic = (truncated.diff() < 0).sum().sum()
        print(f"rank={rank}, nonmonotonic={nonmonotonic}")
