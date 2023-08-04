#!/usr/bin/env python

import numpy as np
import pandas as pd
from pyminimc import plotting

with open(
    "/Users/atumulak/Developer/minimc/benchmarks/continuous_temperature_sensitivity_long.out"
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

prefix = "/Users/atumulak/Developer/minimc/data/tnsl/endfb8/"
plotting.plot_sensitivities(
    mean,
    stddev,
    pd.DataFrame(pd.read_hdf(prefix + "log_offset_beta_S_coarse.hdf5")).shape,
    pd.DataFrame(pd.read_hdf(prefix + "log_offset_beta_CDF_coarse.hdf5"))
    .unstack()
    .shape,
    pd.DataFrame(pd.read_hdf(prefix + "log_offset_beta_E_T_coarse.hdf5"))
    .unstack()
    .shape,
    pd.DataFrame(pd.read_hdf(prefix + "log_alpha_S_coarse.hdf5")).shape,
    pd.DataFrame(pd.read_hdf(prefix + "log_alpha_CDF_coarse.hdf5"))
    .unstack()
    .shape,
    pd.DataFrame(pd.read_hdf(prefix + "log_alpha_beta_T_coarse.hdf5"))
    .unstack()
    .shape,
)
