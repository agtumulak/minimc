#!/usr/bin/env python

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from pyminimc import util, estimator, compression
import pickle

# load reference MCNP energy spectrum
mcnp = estimator.MCNP(
    "/Users/atumulak/Developer/mcnp6-runs/single_zone/single_zone.mctal"
).marginalize("energy")
_, mcnp_densities, _ = mcnp.as_midpoint_density_pairs()


partitions = compression.split(
    pd.DataFrame(pd.read_hdf("data/tsl/endfb8/dfs/beta_full.hdf5")),
    "E",
    [95, 235, 500],
)
truncated_partitions = [
    compression.truncate(partition, rank)
    for partition, rank in zip(partitions, [28, 22, 15, 13])
]
monotonic_partitions = [
    compression.remove_nonmonotonic_cdfs(subset)
    for subset in truncated_partitions
]
with open("trash_bin/coarsened_partitions.pkl", "rb") as f:
    coarsened_partitions = pickle.load(f)
rel_tolerances, rel_l2_errors, worst_indices, rel_peak_errors, cov_sizes_GB = (
    [],
    [],
    [],
    [],
    [],
)
for rel_tolerance in np.linspace(4.6e-4, 1e-3, 55):
    coarsened_partitions = compression.adaptive_coarsen(
        partitions, coarsened_partitions, rel_frobenius_norm_tol=rel_tolerance
    )
    minimc = util.run_minimc_with_tsl(
        "/Users/atumulak/Developer/minimc/build/src/runminimc",
        "/Users/atumulak/Developer/minimc/benchmarks/single_zone.xml",
        coarsened_partitions,
        [
            pd.DataFrame(
                pd.read_hdf(
                    f"/Users/atumulak/Developer/minimc/data/tsl/endfb8/dfs/alpha_{i}.hdf5"
                )
            )
            for i in range(4)
        ],
    ).marginalize("energy")
    # compute errors
    _, minimc_densities, _ = minimc.as_midpoint_density_pairs()
    rel_tolerances.append(rel_tolerance)
    print(f"rel_tolerance: {rel_tolerance}")
    # compute l2 error
    rel_l2_error = np.sum(
        np.subtract(minimc_densities, mcnp_densities[1:])
        * np.diff(minimc._boundaries)
    ) / np.sum(mcnp._counts[1:])
    rel_l2_errors.append(rel_l2_error)
    print(f"rel_l2_error: {rel_l2_error}")
    # compute error in spectrum peak
    i_max = np.argmax(minimc_densities)
    worst_indices.append(i_max)
    print(f"worst index: {i_max}")
    i_max = 6
    rel_peak_error = (
        minimc_densities[i_max] - mcnp_densities[1:][i_max]
    ) / mcnp_densities[1:][6]
    rel_peak_errors.append(rel_peak_error)
    print(f"rel_peak_error: {rel_peak_error}")
    # compute covariance size in GB
    n_params = sum(
        (partition.index.size + 1 + partition.columns.size)
        * np.linalg.matrix_rank(partition)
        for partition in coarsened_partitions
    )
    cov_size_GB = n_params * (n_params + 1) / 2 * 8 / 1000000000
    cov_sizes_GB.append(cov_size_GB)
    print(f"size: {cov_size_GB} GB")

plt.legend()
plt.show()
