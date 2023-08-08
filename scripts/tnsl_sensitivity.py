#!/usr/bin/env python

import copy
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET
from pyminimc import plotting, util
from ipdb import set_trace as st


if __name__ == "__main__":
    # loop over each matrix, perturb it, and compare against sensitivity result
    pd.set_option("display.precision", 8)
    dleak = pd.DataFrame(
        0.0,
        columns=[
            "index",
            "FD (mean)",
            "FD (std. dev.)",
            "DOS (mean)",
            "DOS (std. dev.)",
        ],
        dtype=float,
        index=pd.MultiIndex.from_product([["beta", "alpha"], ["S", "U", "V"]]),
    )
    dleak["index"] = dleak["index"].astype(int)

    with open(
        "/Users/atumulak/Developer/minimc/benchmarks/continuous_temperature_sensitivity.out"
    ) as f:
        line = f.readline()
        while not line.startswith("mean"):
            line = f.readline()
        line = f.readline()
        unperturbed_mean = float(f.readline().split(",")[0])
        while not line.startswith("std dev"):
            line = f.readline()
        line = f.readline()
        unperturbed_stddev = float(f.readline().split(",")[0])
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

    # structure inputs
    prefix = "/Users/atumulak/Developer/minimc/data/tnsl/endfb8/"
    tnsl: util.TNSLType = {}
    offsets = {}
    offset = 0
    for dataset, dataset_prefix, V_matrix_name in [
        ("beta", "log_offset_beta", "E_T"),
        ("alpha", "log_alpha", "beta_T"),
    ]:
        tnsl[dataset] = {}
        offsets[dataset] = {}
        for matrix, matrix_filename in [
            ("S", "S"),
            ("U", "CDF"),
            ("V", V_matrix_name),
        ]:
            file_extension = ".hdf5"
            filepath = f"{prefix}{dataset_prefix}_{matrix_filename}_coarse{file_extension}"
            unperturbed = pd.DataFrame(pd.read_hdf(filepath))
            tnsl[dataset][matrix] = unperturbed
            offsets[dataset][matrix] = offset
            offset += tnsl[dataset][matrix].size
    # most significant element from each matrix
    for dataset in tnsl:
        for matrix in tnsl[dataset]:
            # get flattend indices relevant to current matrix
            begin = offsets[dataset][matrix]
            end = begin + tnsl[dataset][matrix].size
            # search sensitivities to find which parameter has the largest
            # statistical deviation (z score) from zero
            sensitivity_mean = mean[begin:end]
            sensitivity_stddev = stddev[begin:end]
            zscores = sensitivity_mean / sensitivity_stddev
            # create a copy of dataset with most significant parameter
            # perturbed by 5% in the positive direction
            perturbed_tnsl = copy.deepcopy(tnsl)
            param = perturbed_tnsl[dataset][matrix].values.flatten()
            dparam = 0.05 * np.abs(param)
            perturbed_param = param + dparam
            destimate_mean = sensitivity_mean * dparam
            destimate_stddev = sensitivity_stddev * dparam
            # only choose perturbations whose effect is above statistical noise
            # of forward differencing
            is_fd_detectable = np.abs(destimate_mean) > 0.0005
            if is_fd_detectable.sum() == 0:
                print("nothing is forward-differencing detectable")
                st()
            candidate_quality = np.abs(
                np.nan_to_num(np.where(is_fd_detectable, zscores, np.nan))
            )
            best_i = np.argmax(candidate_quality)
            perturbed_tnsl[dataset][matrix].iloc[best_i] = perturbed_param[
                best_i
            ]
            # compute change in leakage using DOS
            row = (dataset, matrix)
            dleak.loc[row, "index"] = best_i
            dleak.loc[row, "DOS (mean)"] = destimate_mean[best_i]
            dleak.loc[row, "DOS (std. dev.)"] = destimate_stddev[best_i]
            # run perturbed dataset and get perturbed leakge and then
            # immediately undo perturbation
            lines = util.run_minimc_tnsl(
                "/Users/atumulak/Developer/minimc/build/src/runminimc",
                "/Users/atumulak/Developer/minimc/benchmarks/continuous_temperature.xml",
                perturbed_tnsl,
            )
            # parse output
            for i, line in enumerate(lines):
                if line.startswith("mean"):
                    perturbed_mean = float(lines[i + 2].split(",")[0])
                    dleak.loc[row, "FD (mean)"] = (
                        perturbed_mean - unperturbed_mean
                    )
                if line.startswith("std dev"):
                    perturbed_stddev = float(lines[i + 2].split(",")[0])
                    dleak.loc[row, "FD (std. dev.)"] = np.sqrt(
                        perturbed_stddev**2 + unperturbed_stddev**2
                    )
            print(dleak)
    # compute Wasserstein metric
    dleak.to_hdf("dleak.hdf5", "pandas")
    st()
    pass
