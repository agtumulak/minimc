#!/usr/bin/env python

import pandas as pd
import numpy as np
from ipdb import set_trace as st


if __name__ == "__main__":
    df = pd.DataFrame(
        pd.read_hdf(
            "/Users/atumulak/Developer/minimc/data/tnsl/endfb8/dfs/scatter_xs.hdf5"
        )
    )
    # convert energy from eV to MeV
    df.index *= 1e-6

    # coarsen energies
    fine_E_grid = df.index
    coarse_E_grid = fine_E_grid[
        np.linspace(0, fine_E_grid.size - 1, 250).astype(int)
    ]
    df = df.loc[coarse_E_grid]

    # get full SVD
    U, S, Vt = np.linalg.svd(df, full_matrices=False)

    # compute relative Frobenius norm between true matrix and low-rank
    # approximation
    rel_frobenius_norm = np.sqrt(np.square(S)[::-1].cumsum()[::-1]) / np.sqrt(
        np.square(S).sum()
    )
    # get first index below specified value
    order = np.where(rel_frobenius_norm < 1e-5)[0][0]
    print(f"order: {order}")

    # construct SVD DataFrames
    U_df = pd.DataFrame(
        {"coefficient": U[:, :order].flatten()},
        index=pd.MultiIndex.from_product(
            [df.index, range(order)], names=["E", "order"]
        ),
    )

    S_df = pd.DataFrame(
        {"coefficient": S[:order]},
        index=pd.MultiIndex.from_product([range(order)], names=["order"]),
    )

    V_df = pd.DataFrame(
        {"coefficient": Vt[:order, :].T.flatten()},
        index=pd.MultiIndex.from_product(
            [df.columns.unique(0), range(order)],
            names=[
                "T",
                "order",
            ],
        ),
    )

    # print worst relative error
    reconstructed = (
        U_df.unstack().values
        @ np.diag(S_df.values.flatten())
        @ V_df.unstack().T.values
    )
    print(
        f"worst relative error: {((reconstructed - df) / df).abs().max().max()}"
    )

    # check format
    U_df.to_hdf(
        "/Users/atumulak/Developer/minimc/data/tnsl/endfb8/scatter_xs_E.hdf5",
        "pandas",
    )
    S_df.to_hdf(
        "/Users/atumulak/Developer/minimc/data/tnsl/endfb8/scatter_xs_S.hdf5",
        "pandas",
    )
    V_df.to_hdf(
        "/Users/atumulak/Developer/minimc/data/tnsl/endfb8/scatter_xs_T.hdf5",
        "pandas",
    )
