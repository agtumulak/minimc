#!/ usr / bin / env python

import numpy as np
import pandas as pd
from pyminimc import plotting

n_partitions = 4

prefix = "~/Developer/minimc/data/tnsl/endfb8/"
partitions = [
    {
        key: pd.DataFrame(
            pd.read_hdf(f"{prefix}alpha_endfb8_{i}_{filename}_coeffs.hdf5")
        )
        for key, filename in (("S", "S"), ("U", "CDF"), ("V", "beta_T"))
    }
    for i in range(n_partitions)
]

with open(
    "/Users/atumulak/Developer/minimc/build/test/broomstick_sensitivity_smallhist.out"
) as f:
    line = f.readline()
    while not line.startswith("leakage::tnsl"):
        line = f.readline()
    while not line.startswith("mean"):
        line = f.readline()
    line = f.readline()
    mean = np.asarray([float(s) for s in f.readline().split(",")[:-1]])
    while not line.startswith("std dev"):
        line = f.readline()
    line = f.readline()
    stddev = np.asarray([float(s) for s in f.readline().split(",")[:-1]])

plotting.plot_sensitivities(mean, stddev, partitions)
