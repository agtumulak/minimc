#!/usr/bin/env python
"""Preprocesses S(a,b,T) data into CDF lookup tables

Based on:
Andrew T. Pavlou, Wei Ji,
On-the-fly sampling of temperature-dependent thermal neutron scattering data for Monte Carlo simulations,
Annals of Nuclear Energy,
Volume 71,
2014,
Pages 411-426,
ISSN 0306-4549,
https://doi.org/10.1016/j.anucene.2014.04.028.
"""

import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import pickle
import re
from collections import OrderedDict
from dask.diagnostics import ProgressBar
from multiprocessing import Pool
from scipy import interpolate
from tqdm import tqdm


# target nuclide mass ratio
A = 0.999167339

# Boltzmann constant in eV / K
k = 8.617333262145e-5

# bound cross section
sigma_free = 4.095600e1 / 2  # 2 hydrogens per H2O
sigma_bound = sigma_free * ((A + 1) / A) ** 2


def parse_file7(mf7_path):
    """
    Parse ENDF File 7 file and return pandas DataFrame.

    Parameters
    ----------
    mf7_path : string
        Path to File 7A
    """

    def to_float(endf_float_string):
        """Convert ENDF-style float string to float."""
        pattern = re.compile(r"\d([+-])")
        return float(pattern.sub(r"E\1", endf_float_string))

    with open(mf7_path, mode="r") as sab_file:
        # skip headers
        while True:
            if sab_file.readline().endswith("1 7  4\n"):
                break
        # skip to relevant part
        for _ in range(4):
            next(sab_file)
        N_beta = int(sab_file.readline().split()[0])
        alphas, betas, Ts, Ss = [], [], [], []
        for _ in range(N_beta):
            # Read first temperature, beta block
            entries = sab_file.readline().split()
            temp, beta = (to_float(x) for x in entries[:2])
            N_temp = int(entries[2]) + 1
            # The first temperature is handled differently: alpha and
            # S(alpha, beta) values are given together.
            N_alpha = int(sab_file.readline().split()[0])
            N_full_rows, remainder = divmod(N_alpha, 3)
            for _ in range(N_full_rows + (remainder != 0)):
                # Everything after column 66 is ignored
                line = sab_file.readline()
                doubles = [
                    to_float(line[start : start + 11])
                    for start in range(0, 66, 11)
                    if not line[start : start + 11].isspace()
                ]
                for alpha, S in zip(doubles[::2], doubles[1::2]):
                    alphas.append(alpha)
                    betas.append(beta)
                    Ts.append(temp)
                    Ss.append(S)
            # The remaining temperatures are handled uniformly
            N_full_rows, remainder = divmod(N_alpha, 6)
            for _ in range(N_temp - 1):
                temp, beta = (
                    to_float(x) for x in sab_file.readline().split()[:2]
                )
                # Subsequent betas use the first beta's alpha grid.
                unique_alphas = (a for a in alphas[:N_alpha])
                for row in range(N_full_rows + (remainder != 0)):
                    line = sab_file.readline()[:66]
                    for S in [
                        to_float(line[start : start + 11])
                        for start in range(0, 66, 11)
                        if not line[start : start + 11].isspace()
                    ]:
                        alphas.append(next(unique_alphas))
                        betas.append(beta)
                        Ts.append(temp)
                        Ss.append(S)
        df = pd.DataFrame.from_dict(
            {"alpha": alphas, "beta": betas, "T": Ts, "S": Ss}
        )
        # rescale alpha and beta due to LAT flag in LEAPR
        df["alpha"] = df["alpha"] * 293.6 / df["T"]
        df["beta"] = df["beta"] * 293.6 / df["T"]
        df = df.set_index(["beta", "alpha", "T"])
        return df


def lin_log_cum_trapz(s):
    """
    Integrates a Series with trapezoid rule using linear in index and
    logarithmic in value interpolation

    Parameters
    ----------
    s : pd.Series
       Series to integrate
    """
    x, y = s.index, s.values
    numerators = y[1:] - y[:-1]
    denominators = np.log(y[1:]) - np.log(y[:-1])
    # https://stackoverflow.com/a/37977222/5101335
    ratios = np.divide(
        numerators, denominators, out=numerators, where=denominators != 0
    )
    ratios = np.nan_to_num(ratios)
    return pd.Series(
        np.concatenate(([0], np.cumsum(ratios * (x[1:] - x[:-1])))),
        index=s.index,
    )


def process_E_T(args):
    """
    Generates PDFs in beta, conditional PDFs in alpha given beta, and CDFs
    thereof

    Parameters
    ----------
    args : tuple
        (S(a,b,T) DataFrame, Incident energy in eV, Temperature in K)
    """
    sab_df, E, T = args
    T_betas = np.array(sorted(sab_df.loc[:, :, T].index.unique("beta")))
    # valid beta values
    min_beta = -E / (k * T)
    max_beta = 20
    min_beta_index = np.searchsorted(
        T_betas, -min_beta
    )  # as range end, will never include min_beta
    max_beta_index = np.searchsorted(
        T_betas, max_beta
    )  # as range end, will never include max_beta
    betas = np.concatenate(
        (
            [min_beta],
            -np.flip(T_betas[1:min_beta_index]),
            T_betas[:max_beta_index],
            [max_beta],
        )
    )
    alpha_cdfs = {}
    alpha_pdfs = {}
    beta_pdf = pd.Series(0, index=betas)  # pdf is zero at min_beta and max_beta
    for beta in betas[1:-1]:
        # energy-independent alpha distribution
        S_values = sab_df.xs(
            (beta if beta > 0 else -beta, T), level=("beta", "T")
        )["S"]
        # energy-dependent alpha distribution, valid alpha values
        min_alpha = np.square(np.sqrt(E) - np.sqrt(E + beta * k * T)) / (
            A * k * T
        )
        max_alpha = np.square(np.sqrt(E) + np.sqrt(E + beta * k * T)) / (
            A * k * T
        )
        S_at_min_alpha, S_at_max_alpha = np.interp(
            [min_alpha, max_alpha], S_values.index, S_values
        )
        min_alpha_index = np.searchsorted(
            S_values.index, min_alpha, side="right"
        )  # will never include min_alpha
        max_alpha_index = np.searchsorted(
            S_values.index, max_alpha
        )  # as range end, will never include max_alpha
        S_values = pd.concat(
            (
                pd.Series({min_alpha: S_at_min_alpha}),
                S_values.iloc[min_alpha_index:max_alpha_index],
                pd.Series({max_alpha: S_at_max_alpha}),
            )
        )
        alpha_cdf = lin_log_cum_trapz(S_values)
        alpha_integral = alpha_cdf.iloc[-1]
        alpha_cdfs[beta] = alpha_cdf / alpha_integral
        alpha_pdfs[beta] = (S_values / alpha_integral).rename(beta)
        beta_pdf.loc[beta] = np.exp(-beta / 2.0) * alpha_integral
    # convert pdf to cdf
    beta_cdf = lin_log_cum_trapz(beta_pdf)
    beta_integral = beta_cdf.iloc[-1]
    beta_cdf /= beta_integral
    beta_pdf /= beta_integral
    total_inelastic_xs = sigma_bound * A * k * T / (4 * E) * beta_integral
    return beta_pdf, alpha_pdfs, beta_cdf, alpha_cdfs, total_inelastic_xs


def process_b_T(sab_s, max_alpha):
    """
    Generates conditional CDF in alpha given beta and temperature. Returns None
    if there is zero probability of the given beta being sampled.

    Parameters
    ----------
    sab_s : pd.Series
        A pd.Series containing a MultiIndex. The MultiIndex levels are
        temperature `T`, `beta`, and `alpha`. There is only a single value of
        `T` and `beta` while multiple `alpha` values must be present. The
        values are corresponding value of S(a,b,T).
    max_alpha : double
        largest alpha in the dataset
    """
    sab_s.index = sab_s.index.droplevel(["T", "beta"])
    # set endpoints to zero
    sab_s.loc[0] = 0
    sab_s.loc[max_alpha] = 0
    sab_s = sab_s.sort_index()
    E_independent_alpha_cdf = lin_log_cum_trapz(sab_s)
    E_independent_alpha_integral = E_independent_alpha_cdf.iloc[-1]
    # if integral is zero, return nothing
    if E_independent_alpha_integral == 0:
        E_independent_alpha_cdf = None
    else:
        E_independent_alpha_cdf /= E_independent_alpha_integral
    return E_independent_alpha_cdf


def beta_fitting_function(T, c0, c1, c2, c3, c4, c5):
    return (
        c0
        + c1 / T ** (1 / 2)
        + c2 / T
        + c3 / T ** (3 / 2)
        + c4 / T**2
        + c5 / T ** (5 / 2)
    )


def alpha_fitting_function(T, c0, c1, c2):
    return c0 + c1 / T + c2 / T**2


def get_pdf_direct(sab_df, E, T, label=None):
    """
    Creates bivariate PDF in alpha and beta from given S(a,b,T) DataFrame

    Alpha grid is unionized across all possible beta values. Linear
    interpolation is used in alpha for missing values. Zero values are used
    outside valid alpha ranges.

    Parameters
    ----------
    sab_df : pd.DataFrame
        S(a,b,T) DataFrame
    E : float
        Incident energy in eV
    T : float
        Temperature in K
    """
    label = label if label else "direct"
    beta_pdf, alpha_pdfs, _, _, _ = process_E_T((sab_df, E, T))
    for beta, p_beta in beta_pdf.iloc[1:-1].iteritems():
        alpha_pdfs[beta] *= p_beta
    return (
        pd.concat(alpha_pdfs.values(), names=alpha_pdfs.keys(), axis="columns")
        .interpolate(method="index", axis="index", limit_area="inside")
        .fillna(0)
        .stack()
        .rename_axis(["alpha", "beta"])
        .rename(label)
    )


def get_pdf_pod(
    beta_T_path,
    beta_S_path,
    beta_E_CDF_path,
    alpha_T_path,
    alpha_S_path,
    alpha_beta_CDF_path,
    E,
    T,
    max_alpha,
    label=None,
):
    """
    Reconstructs bivariate PDF in alpha and beta from proper orthogonal
    decomposition coefficients for alpha and beta

    Parameters
    ---------
    beta_T_path : string
        Path to HDF5 file containing temperature dependent coefficients for beta
    beta_S_path : string
        Path to HDF5 file containing singular values for beta
    beta_E_CDF : string
        Path to HDF5 file containing energy and CDF dependent coefficients for
        beta
    alpha_T_path : string
        Path to HDF5 file containing temperature dependent coefficients for alpha
    alpha_S_path : string
        Path to HDF5 file containing singular values for alpha
    alpha_E_CDF : string
        Path to HDF5 file containing beta and CDF dependent coefficients for
        alpha
    E : float
        Incident energy in eV
    T : float
        Temperature in K
    max_alpha : float
        Largest alpha in the dataset
    """
    label = label if label else "direct (pod reconstructed)"
    # Load data
    beta_T = pd.read_hdf(beta_T_path)
    beta_S = pd.read_hdf(beta_S_path)
    beta_E_CDF = pd.read_hdf(beta_E_CDF_path)
    alpha_T = pd.read_hdf(alpha_T_path)
    alpha_S = pd.read_hdf(alpha_S_path)
    alpha_beta_CDF = pd.read_hdf(alpha_beta_CDF_path)
    # evaluate beta at nearest E and nearest T
    Ts = beta_T.index.unique("T")
    nearest_T = Ts[np.argmin(np.abs(Ts - T))]
    Es = beta_E_CDF.index.unique("E")
    nearest_E = Es[np.argmin(np.abs(Es - E * 1e-6))]
    # Reconstruct beta CDF
    beta_cdf = pd.Series(
        (
            beta_T.loc[nearest_T].T.values
            @ np.diag(beta_S.values.flatten())
            @ beta_E_CDF.loc[nearest_E].unstack().T.values
        ).flatten(),
        index=beta_E_CDF.index.unique("CDF"),
    )
    # compute PDF
    beta_pdf = pd.Series(
        (beta_cdf.index[1:] - beta_cdf.index[:-1])
        / (beta_cdf.values[1:] - beta_cdf.values[:-1]),
        index=beta_cdf.values[:-1],
    )
    # evaluate alpha at (possibly different) nearest T and nearest beta
    Ts = alpha_T.index.unique("T")
    nearest_T = Ts[np.argmin(np.abs(Ts - T))]
    # Reconstruct alpha CDF
    alpha_cdf = pd.Series(
        (
            alpha_T.loc[nearest_T].T.values
            @ np.diag(alpha_S.values.flatten())
            @ alpha_beta_CDF.unstack().T.values
        ).flatten(),
        index=alpha_beta_CDF.unstack().index,
    ).unstack("beta")
    # append alpha CDFs for negative beta values
    alpha_betas = alpha_cdf.columns
    # find largest beta in alpha_betas which is strictly less than E / (k * T)
    # we assume beta = 0 exists so result of searchsorted is >= 1
    min_beta = alpha_betas[np.searchsorted(alpha_betas, -beta_cdf.iloc[0]) - 1]
    neg_b_alpha_cdf = (
        alpha_cdf.loc[:, alpha_betas[1] : min_beta]  # don't include beta = 0
        .rename(columns=lambda x: -x)  # make beta labels negative
        .sort_index(axis="columns")
    )
    # find largest beta in alpha_betas which is strictly less than 20
    max_beta = alpha_betas[np.searchsorted(alpha_betas, 20) - 1]
    alpha_cdf = pd.concat(
        (neg_b_alpha_cdf, alpha_cdf.loc[:max_beta]), axis="columns"
    )
    # add endpoints
    alpha_cdf.loc[0, :] = 0
    alpha_cdf.loc[1, :] = max_alpha
    alpha_cdf = alpha_cdf.sort_index()
    # choose common alpha_grid
    alphas = np.linspace(0, max_alpha, 10000)
    # choose subset of alpha values
    def get_joint_probability(s):
        """
        Multiplies conditional probability in alpha with probability in beta
        """
        nonlocal alphas
        nonlocal E
        beta = s.name
        # skip values of beta which won't be reached
        if E + beta * k * T < 0:
            return pd.Series(0, index=alphas)
        # choose correct subset of CDF values within min_alpha and max_alpha
        s = pd.Series(s.index, index=pd.Index(s, name="alpha"))
        # insert values for min_alpha and max_apha
        min_alpha = np.square(np.sqrt(E) - np.sqrt(E + beta * k * T)) / (
            A * k * T
        )
        max_alpha = np.square(np.sqrt(E) + np.sqrt(E + beta * k * T)) / (
            A * k * T
        )
        s.loc[min_alpha] = np.nan  # interpolate value later
        s.loc[max_alpha] = np.nan  # interpolate value later

        s = s.sort_index().interpolate(method="index").loc[min_alpha:max_alpha]
        # rescale CDF to be 0 at min_alpha and 1 at max_alpha
        s = (s - s.min()) / (s.max() - s.min())
        alpha_pdf = pd.Series(
            (s.values[1:] - s.values[:-1]) / (s.index[1:] - s.index[:-1]),
            index=s.index[:-1],
        )
        alpha_pdf.loc[min_alpha] = 0
        alpha_pdf.loc[max_alpha] = 0
        alpha_pdf = alpha_pdf.sort_index()
        # use common alpha grid
        new_alphas = set(alphas).difference(set(alpha_pdf.index))
        alpha_pdf = (
            pd.concat([alpha_pdf, pd.Series(np.nan, index=new_alphas)])
            .sort_index()
            .interpolate(method="index", limit_area="inside")
            .fillna(0)
            .loc[alphas]
        )
        # interpolate value of beta pdf
        nonlocal beta_pdf
        beta_pdf_value = np.interp(
            beta, beta_pdf.index, beta_pdf, left=0, right=0
        )
        return alpha_pdf * beta_pdf_value

    bivariate_pdf = alpha_cdf.apply(get_joint_probability).stack()
    bivariate_pdf.index.names = ["alpha", "beta"]
    return bivariate_pdf.rename(label)


def get_pdf_runsab(counts_path, *bounds_paths):
    """
    Creates bivariate PDF in alpha and beta from a list of counts and N paths
    to axis boundaries

    Parameters
    ----------
    counts_path : string
        Path to file containing counts separated by newlines. Counts are
        ordered with the last element in `bounds_paths` changing the fastest,
        and the first element in `bounds_paths` changing slowest.
    bounds_paths : sequence of strings
        paths to files containing bin boundaries.

    Returns
    -------
    pd.Series with a MultiIndex for each axis
    """
    with open(counts_path) as f:
        counts = np.array([int(l.strip()) for l in f.readlines()])
        counts = counts / counts.sum()
    widths, bin_edges = [], []
    for bounds_path in bounds_paths:
        with open(bounds_path) as f:
            bounds = np.array([float(x) for x in f.readlines()])
            widths.append(bounds[1:] - bounds[:-1])
            bin_edges.append(bounds[:-1])
    bin_areas = np.einsum("i,j->ij", *widths).reshape(-1)
    density = counts / bin_areas
    return pd.Series(
        density,
        index=pd.MultiIndex.from_product(bin_edges, names=["alpha", "beta"]),
        name="minimc",
    )


def get_pdf_mcnp(mctal_path, E, T, label=None):
    """
    Creates bivariate PDF in alpha and beta from an MCNP mctal file

    Parameters
    ----------
    mctal_path : string
        Path to mctal file
    E : float
        Incident energy in eV
    T : float
        Temperature in K
    """
    label = label if label else "mcnp"
    cosine_bounds = [
        -1.0,
    ]
    energy_bounds = [
        0.0,
    ]
    counts = []
    with open(mctal_path) as f:
        line = ""
        # skip to cosine bins
        while not line.startswith("c"):
            line = f.readline()
        line = f.readline()
        # parse cosine boundaries until energy boundaries section is reached
        while not line.startswith("e"):
            cosine_bounds.extend(float(x) for x in line.split())
            line = f.readline()
        line = f.readline()
        # parse energy boundaries until time boundaries section is reached
        while not line.startswith("t"):
            energy_bounds.extend(float(x) * 1e6 for x in line.split())
            line = f.readline()
        # skip to values
        while not line.startswith("vals"):
            line = f.readline()
        line = f.readline()
        # parse values until the tfc section is reached
        while not line.startswith("tfc"):
            counts.extend(float(x) for x in line.split()[::2])
            line = f.readline()
    # compute densities in cosine, eV space
    cosine_bounds = np.array(cosine_bounds)
    energy_bounds = np.array(energy_bounds)
    counts = np.array(counts)
    cosine_widths = cosine_bounds[1:] - cosine_bounds[:-1]
    energy_widths = energy_bounds[1:] - energy_bounds[:-1]
    bin_areas = np.einsum("i,j->ij", cosine_widths, energy_widths).reshape(-1)
    density = counts / bin_areas
    s = pd.Series(
        density,
        index=pd.MultiIndex.from_product(
            [cosine_bounds[:-1], energy_bounds[:-1]], names=["mu", "E"]
        ),
        name="mcnp",
    )

    def to_alpha_beta(s):
        beta = s.name
        out_E = E + beta * k * T
        mus = s.index.get_level_values("mu")
        alphas = (E + out_E - 2 * mus * np.sqrt(E * out_E)) / (A * k * T)
        return pd.Series(s.values, index=alphas).rename_axis(
            index={"mu": "alpha"}
        )

    # multiply by jacobian
    s = (
        s
        * 0.5
        * A
        * (k * T) ** 2
        / (np.sqrt(s.index.get_level_values("E") * E))
    )
    s = (
        s.rename(index=lambda out_E: (out_E - E) / (k * T), level="E")
        .rename_axis(index={"E": "beta"})
        .groupby("beta")
        .apply(to_alpha_beta)
    )
    return (
        s[~s.index.duplicated()]
        .unstack("beta")
        .interpolate(method="index", axis="index", limit_area="inside")
        .fillna(0)
        .stack()
        .rename_axis(["alpha", "beta"])
        .rename(label)
    )


def get_pdf_minimc(minimc_path, E, T, label=None):
    """
    Creates bivariate PDF in alpha and beta from minimc output
    """
    label = label if label else "minimc"
    with open(minimc_path) as f:
        line = ""
        # skip to cosine bins
        while not line.startswith("cosine"):
            line = f.readline()
        line = f.readline()
        cosine_bounds = [float(x) for x in f.readline().split(",")[:-1]]
        # skip to energy bins
        while not line.startswith("energy"):
            line = f.readline()
        line = f.readline()
        energy_bounds = [float(x) * 1e6 for x in f.readline().split(",")[:-1]]
        # skip to values
        while not line.startswith("mean"):
            line = f.readline()
        line = f.readline()
        counts = [float(x) for x in f.readline().split(",")[:-1]]
    # compute densities in cosine, eV space
    cosine_bounds = np.concatenate(([-np.inf], cosine_bounds, [np.inf]))
    energy_bounds = np.concatenate(([-np.inf], energy_bounds, [np.inf]))
    counts = np.array(counts)
    cosine_widths = cosine_bounds[1:] - cosine_bounds[:-1]
    energy_widths = energy_bounds[1:] - energy_bounds[:-1]
    bin_areas = np.einsum("i,j->ij", cosine_widths, energy_widths).reshape(-1)
    # compute densities and remove infinitely sized bins
    density = (
        (counts / bin_areas)
        .reshape(cosine_widths.size, energy_widths.size)[1:-1, 1:-1]
        .reshape(-1)
    )
    # create pd.Series
    cosine_midpoints = (cosine_bounds[2:-1] + cosine_bounds[1:-2]) / 2
    energy_midpoints = (energy_bounds[2:-1] + energy_bounds[1:-2]) / 2
    s = pd.Series(
        density,
        index=pd.MultiIndex.from_product(
            [cosine_midpoints, energy_midpoints], names=["mu", "E"]
        ),
        name="minimc",
    )
    # convert to alpha, beta space; #TODO: this is identical to mcnp version
    def to_alpha_beta(s):
        beta = s.name
        out_E = E + beta * k * T
        mus = s.index.get_level_values("mu")
        alphas = (E + out_E - 2 * mus * np.sqrt(E * out_E)) / (A * k * T)
        return pd.Series(s.values, index=alphas).rename_axis(
            index={"mu": "alpha"}
        )

    # multiply by jacobian
    s = (
        s
        * 0.5
        * A
        * (k * T) ** 2
        / (np.sqrt(s.index.get_level_values("E") * E))
    )
    s = (
        s.rename(index=lambda out_E: (out_E - E) / (k * T), level="E")
        .rename_axis(index={"E": "beta"})
        .groupby("beta")
        .apply(to_alpha_beta)
    )
    return (
        s[~s.index.duplicated()]
        .unstack("beta")
        .interpolate(method="index", axis="index", limit_area="inside")
        .fillna(0)
        .stack()
        .rename_axis(["alpha", "beta"])
        .rename(label)
    )


def plot_pdf(s):
    """
    Plots bivariate PDF of in alpha and beta

    Parameters
    ----------
    s : pd.Series
        Bivariate PDF in alpha and beta to plot. MultiIndex must be beta
        followed by alpha.
    """
    plt.contourf(
        s.index.unique("beta"),
        s.index.unique("alpha"),
        np.log(s).unstack(),
        levels=100,
    )
    plt.xlabel(r"$\beta$")
    plt.ylabel(r"$\alpha$")
    plt.title(r"$\log p_{\alpha, \beta} (\alpha, \beta)$")
    plt.colorbar()
    plt.show()


def compare_bivariate_pdf(title, *series):
    """
    Compares PDFs in alpha and beta from multiple series.

    Parameters
    ----------
    title : string
        Plot title
    s1, s2, ... : sequence of PDFs in alpha and beta
    """
    nrows, ncols = 2, 2
    min_log_density = np.log(min(s[s > 0].min() for s in series))
    min_beta = max(s.index.unique("beta").min() for s in series)
    max_beta = min(s.index.unique("beta").max() for s in series)
    min_alpha = max(s.index.unique("alpha").min() for s in series)
    max_alpha = min(s.index.unique("alpha").max() for s in series)
    f, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        sharex="col",
        sharey="row",
        constrained_layout=True,
    )
    for i, row in enumerate(axs):
        for j, ax in enumerate(row):
            k = ncols * i + j
            if k >= len(series):
                break
            s = series[k]
            cm = ax.contourf(
                s.index.unique("beta"),
                s.index.unique("alpha"),
                np.log(s).unstack(),
                levels=np.linspace(min_log_density, 0, 100),
            )
            if j == 0:
                ax.set_ylabel(r"$\alpha$")
            if i == 1:
                ax.set_xlabel(r"$\beta$")
            ax.set_xlim(min_beta, max_beta)
            ax.set_ylim(min_alpha, max_alpha)
            ax.set_title(s.name)
    f.colorbar(cm, ax=axs, location="right")
    plt.show()


def marginalize(s, axis="beta"):
    """
    Converts a bivariate PDF into a univariate PDF in given axis

    Parameters
    ----------
    s : pd.Series
        Bivariate PDF in alpha and beta to plot. MultiIndex must be beta
        followed by alpha.
    axis : {'alpha', 'beta'}, optional
        The axis of the resulting univariate PDF
    """
    return s.groupby(axis).apply(
        lambda s: np.trapz(s.values, s.index.droplevel(axis))
    )


def compare_univariate_pdf(title, *series, axis="beta"):
    """
    Compares PDFs in beta from multiple series.

    Parameters
    ----------
    title : string
        Plot title
    s1, s2, ... : sequence of PDFs in beta
        Beta PDFs to plot
    axis : {'alpha', 'beta'}, optional
        The axis corresponding to the abscissa
    """
    for s in series:
        s.plot(label=s.name)
    plt.xlabel(rf"$\{axis}$")
    plt.ylabel(rf"$p_{{\{axis}}} (\{axis})$")
    plt.xlim(
        max(s.index.min() for s in series), min(s.index.max() for s in series)
    )
    plt.legend()
    plt.title(title)
    plt.show()


def beta_functional_expansion(
    sab_df, E_min=1e-5, E_max=4.0, n_Es=1000, n_cdfs=1000, order=None
):
    """
    Computes the CDF in beta at various incident energies and temperatures,
    then performs a functional expansion in temperature at various incident
    energies and CDF values.

    Parameters
    ----------
    sab_df : pd.DataFrame
        S(a,b,T) DataFrame
    E_min : float, optional
        Minimum incident energy in eV. Will be included in final energy grid.
    E_max : float, optional
        Maximum incident energy in eV. Will not be included in final energy
        grid.
    n_Es : int, optional
        Approximate number of incident energies (equally spaced in lethargy)
    n_cdfs : int, optional
        Number of CDF values to use
    order : int, optional
        Expansion order for proper orthogonal decomposition. Setting to None
        will return the full expansion.

    Returns
    -------
    pd.DataFrame
        Expansion coefficients in temperature for beta given a given set of
        incident energies and CDF values.

    Todo
    ----
    Use DataFrames for intermediate steps instead of numpy arrays. This is how
    alpha_functional_expansion() does it.
    """
    # Populate `beta_cdfs`, a 2D array where the first index corresponds to
    # incident energy and the second index corresponds to temperature. Each
    # element of this array is a cumulative distribution function for beta.
    df_Ts = np.array(sorted(sab_df.index.unique("T")))
    # equally spaced-lethargy intervals, do not include zero lethargy
    assert E_min > 0
    lethargies = np.linspace(0, np.log(E_max / E_min), num=n_Es + 1)[:0:-1]
    Es = E_max * np.exp(-lethargies)
    with Pool(processes=8) as pool:
        # incident energy, temperature pairs
        E_T_values = np.array([(sab_df, E, T) for E in Es for T in df_Ts])
        results = np.array(
            [
                [x[2], x[4]]
                for x in tqdm(
                    pool.imap(func=process_E_T, iterable=E_T_values),
                    total=len(E_T_values),
                )
            ],
            dtype=object,
        ).reshape(len(Es), len(df_Ts), -1)
        beta_cdfs = results[:, :, 0]
        inelastic_xs = pd.DataFrame(results[:, :, 1], index=Es, columns=df_Ts)
    F = np.linspace(0, 1, n_cdfs)
    beta_df = pd.DataFrame(
        np.nan,
        index=pd.Index(F, name="CDF"),
        columns=pd.MultiIndex.from_product((Es, df_Ts), names=("E", "T")),
    )
    # a good beta grid for the alpha functional expansion
    unique_betas = np.unique(np.abs(beta_df))
    N_betas = 1000
    selected_betas = unique_betas[
        np.round(np.linspace(0, unique_betas.size - 1, N_betas)).astype(int)
    ]
    # interpolate to fill in selected CDF values
    for E, x_E in zip(Es, beta_cdfs):
        for T, beta_cdf in zip(df_Ts, x_E):
            beta_df.loc[:, (E, T)] = np.interp(F, beta_cdf, beta_cdf.index)
    # perform proper orthogonal decomposition
    beta_df_pod_form = beta_df
    U, S, Vt = np.linalg.svd(beta_df_pod_form, full_matrices=False)
    order = S.size if order is None else order
    U_df = pd.DataFrame(
        {"coefficient": U[:, :order].flatten()},
        index=pd.MultiIndex.from_product(
            [F, range(order)], names=["CDF", "order"]
        ),
    )
    S_df = pd.DataFrame(
        {"coefficient": S[:order]},
        index=pd.MultiIndex.from_product([range(order)], names=["order"]),
    )
    V_df = pd.DataFrame(
        {"coefficient": Vt.T[:, :order].flatten()},
        index=pd.MultiIndex.from_product(
            [Es, df_Ts, range(order)], names=["E", "T", "order"]
        ),
    )
    # set energy units to MeV
    V_df.index = V_df.index.set_levels(
        V_df.index.unique(level="E") * 1e-6, level="E"
    )
    # reconstruct
    beta_df_reconstructed = pd.DataFrame(
        U[:, :order] @ np.diag(S[:order]) @ Vt[:order, :],
        index=beta_df_pod_form.index,
        columns=beta_df_pod_form.columns,
    )
    # check that CDFS are monotonic for certain T values
    print(
        f"RMSE: {np.sqrt(((beta_df_reconstructed - beta_df_pod_form)**2).mean().mean())}"
    )
    is_monotonic = beta_df_reconstructed.apply(lambda s: s.is_monotonic)
    print(
        f"Of {Es.size} incident energies and {df_Ts.size} target "
        f"temperatures, {is_monotonic.sum()} of {is_monotonic.size} "
        f"({is_monotonic.sum() / is_monotonic.size * 100}%) have "
        f"monotonic beta as a function of CDF"
    )
    if not is_monotonic.all():
        print("The following CDFs are not monotonic:")
        print(beta_df_reconstructed.loc[:, ~is_monotonic])
    return U_df, S_df, V_df


def alpha_functional_expansion(sab_df, selected_betas, n_cdfs=1000, order=None):
    """
    Computes the energy-independent conditional CDF in alpha given beta at
    various beta values and temperatures, then performs a functional expansion
    in CDF at various beta and temperature values.

    Parameters
    ----------
    sab_df : pd.DataFrame
        S(a,b,T) DataFrame
    selected_betas : np.ndarray
        Beta values to use
    n_cdfs : int, optional
        Number of CDF values to use
    order : int, optional
        Expansion order for proper orthogonal decomposition. Setting to None
        will return the full expansion.
    """
    df_Ts = np.array(sorted(sab_df.index.unique("T")))
    selected_betas = set(selected_betas)

    def common_beta_grid(group):
        """
        Modifies beta values at this temperature to conform with common beta
        grid.
        """
        nonlocal sab_df
        group.index = group.index.droplevel("T")
        # add betas from selected_betas which are not already in group
        new_betas = selected_betas.difference(group.index.unique("beta"))
        new_df = pd.DataFrame(
            np.nan,
            columns=group.index.unique("alpha"),
            index=pd.Index(new_betas, name="beta"),
        )
        combined_df = pd.concat(
            [group.unstack("alpha")["S"], new_df]
        ).sort_index()
        # S(a,b) above maximum beta at this temperature is zero
        combined_df.loc[
            combined_df.index > group.index.unique("beta").max()
        ] = 0
        # interpolate linearly in beta and linearly in ln(S) (interpolation
        # scheme 4 in ENDF)
        return (
            np.exp(
                np.log(combined_df).interpolate(
                    method="index", axis="index", limit_area="inside"
                )
            )
            .loc[list(selected_betas)]
            .stack()
        )

    common_beta_sab_df = sab_df.groupby("T").apply(common_beta_grid)
    # alpha grids do not have to match across T-beta pairs, but they do have to
    # have the same minimum and maximum values
    print("computing alpha CDFs...")
    largest_alpha = sab_df.index.unique("alpha").max()
    print(f"largest alpha: {largest_alpha}")
    # choose number of CDF points we want to use; don't include 0 or 1
    F = np.linspace(0, 1, n_cdfs + 2)[1:-1]
    # Some T-beta values will have all-zeros and result in missing entries for
    # such T-beta values.
    alpha_df = (
        common_beta_sab_df.groupby(["T", "beta"])
        .apply(process_b_T, largest_alpha)
        .groupby(["T", "beta"])
        .apply(
            lambda s: pd.Series(
                np.interp(F, s, s.index.get_level_values("alpha")),
                index=pd.Index(F, name="CDF"),
            )
        )
    )
    # construct matrix to decompose with POD
    alpha_df_pod_form = alpha_df.unstack(["beta", "T"])
    # add back T-beta pairs which didn't have a well-defined alpha CDF
    undefined_cdfs = common_beta_sab_df.unstack("alpha").index.difference(
        alpha_df.unstack("CDF").index
    )
    alpha_df_pod_form[undefined_cdfs] = np.nan
    alpha_df_pod_form = alpha_df_pod_form.sort_index(axis="columns")

    # impute NaN entries: https://stackoverflow.com/a/35611142/5101335
    nan_entries = alpha_df_pod_form.isna()
    alpha_df_pod_form = alpha_df_pod_form.T.fillna(
        value=alpha_df_pod_form.T.mean()
    ).T
    # repeat multiple times until convergence
    converged = False
    prev_trace_norm = np.nan
    while not converged:
        U, S, Vt = np.linalg.svd(alpha_df_pod_form, full_matrices=False)
        order = S.size if order is None else order
        A = U[:, :order] @ np.diag(S[:order]) @ Vt[:order, :]
        alpha_df_pod_form = alpha_df_pod_form.where(~nan_entries, A)
        trace_norm = S.sum()
        abs_rel_diff = np.abs((trace_norm - prev_trace_norm) / prev_trace_norm)
        print(f"trace norm abs rel diff: {abs_rel_diff}", end=" " * 10 + "\r")
        if abs_rel_diff < 1e-8:
            break
        else:
            prev_trace_norm = trace_norm
    U_df = pd.DataFrame(
        {"coefficient": U[:, :order].flatten()},
        index=pd.MultiIndex.from_product(
            [F, range(order)], names=["CDF", "order"]
        ),
    )
    S_df = pd.DataFrame(
        {"coefficient": S[:order]},
        index=pd.MultiIndex.from_product([range(order)], names=["order"]),
    )
    V_df = pd.DataFrame(
        {"coefficient": Vt.T[:, :order].flatten()},
        index=pd.MultiIndex.from_product(
            [alpha_df.index.unique("beta"), df_Ts, range(order)],
            names=["beta", "T", "order"],
        ),
    )
    # reconstruct
    alpha_df_reconstructed = pd.DataFrame(
        U[:, :order] @ np.diag(S[:order]) @ Vt[:order, :],
        index=alpha_df_pod_form.index,
        columns=alpha_df_pod_form.columns,
    )
    print(
        f"RMSE: {np.sqrt(((alpha_df_reconstructed - alpha_df_pod_form)**2).mean().mean())}"
    )
    # check that CDFS are monotonic for certain T values
    is_monotonic = alpha_df_reconstructed.apply(lambda s: s.is_monotonic)
    print(
        f"{is_monotonic.sum()} of {is_monotonic.size} "
        f"({is_monotonic.sum() / is_monotonic.size * 100}%) have "
        f"monotonic alpha as a function of CDF"
    )
    if not is_monotonic.all():
        print("The following CDFs are not monotonic:")
        print(alpha_df_reconstructed.loc[:, ~is_monotonic])
    return U_df, S_df, V_df


def adaptive_coarsen(
    full_df, rel_err_tol=1e-2, abs_err_tol=1e-12, graphics=False
):
    # removed axes from least important to most important
    removed = []
    # initialize data which will not be updated by coarsening loop
    fine_axes = [x.values for x in full_df.index.levels]
    fine_array = full_df.values.reshape(full_df.index.levshape)
    # initialize data which will be updated by coarsening loop
    coarse_axes = fine_axes.copy()
    coarse_array = fine_array.copy()
    fine_interped = fine_array.copy()
    residuals = fine_interped - fine_array
    local_sqr_errs = (residuals * residuals) / residuals.size
    global_err = 0.0
    fine_idx = [np.arange(axis.size) for axis in fine_axes]
    # prepare figure
    if graphics:
        figure, ax = plt.subplots()
        fig = plt.figure()
        ax = fig.gca()
        im = ax.imshow(fine_array[:, :, 8])
        ax.set_xlabel("CDF Index")
        ax.set_ylabel("Energy Index")
        ax.set_title(r"$\beta$ at 450K")
        figure.show()
    while True:
        # identify least important gridpoint
        selected_axis, selected_idx = None, None
        selected_global_err = np.inf
        for axis in range(coarse_array.ndim):
            fine_axis = fine_axes[axis]
            axis_fine_idx = fine_idx[axis]
            # get contribution of each axis gridpoint to total error
            axis_sqr_err = local_sqr_errs.sum(
                axis=tuple(i for i in range(coarse_array.ndim) if i != axis)
            )
            resize_shapes = [
                1 if a != axis else -1 for a in range(coarse_array.ndim)
            ]
            for coarse_idx in range(1, coarse_array.shape[axis] - 1):
                # get the fine axis indices that will be used as the left/right
                # boundaries for the interpolated points
                left_fine_idx = axis_fine_idx[coarse_idx - 1]
                right_fine_idx = axis_fine_idx[coarse_idx + 1]
                # get fine axis indices for new points that have to be
                # interpolated
                interped_fine_idx = np.arange(left_fine_idx + 1, right_fine_idx)
                # compute the interpolation factor for axis points that will be
                # interpolated
                interp_factor = (
                    (fine_axis[interped_fine_idx] - fine_axis[left_fine_idx])
                    / (fine_axis[right_fine_idx] - fine_axis[left_fine_idx])
                ).reshape(resize_shapes)
                # get the previous iteration's interpolated values that will be
                # used as the left/right faces for the interpolated points
                left_face = fine_interped.take([left_fine_idx], axis=axis)
                right_face = fine_interped.take([right_fine_idx], axis=axis)
                # interpolate fine points
                interped_fine = left_face + interp_factor * (
                    right_face - left_face
                )
                # get error contribution from points which are unaffected by
                # removing gridpoint at coarse_idx
                unaffected_err_sum = (
                    axis_sqr_err[: left_fine_idx + 1].sum()
                    + axis_sqr_err[right_fine_idx:].sum()
                )
                # get absolute relative error of interpolated points
                true_values = fine_array.take(interped_fine_idx, axis=axis)
                residuals = interped_fine - true_values
                affected_err = (residuals * residuals) / fine_array.size
                global_err = np.sqrt(unaffected_err_sum + affected_err.sum())
                # # print diagnostics
                # print(
                #         f"axis: {axis}, coarse_idx: {coarse_idx}, "
                #         f"max_new_local_rel_err: {new_local_rel_err.max()}, "
                #         f"global_err: {global_err}")
                abs_rel_errs = np.abs(residuals / true_values)
                abs_errs = np.abs(residuals)
                # automatically reject this coarse gridpoint if it makes any
                # local error too large
                if np.any(
                    (abs_errs > abs_err_tol) & (abs_rel_errs > rel_err_tol)
                ):
                    continue
                # if removing this gridpoint doesn't improve on the candidate
                # gridpoint, skip it
                elif global_err > selected_global_err:
                    continue
                # update
                else:
                    selected_axis = axis
                    selected_idx = coarse_idx
                    selected_global_err = global_err
                    selected_interped = interped_fine
                    selected_interped_fine_idx = interped_fine_idx
        if selected_axis != None:
            # update figure
            if graphics:
                if selected_axis == 0:
                    im.axes.axhline(
                        fine_idx[selected_axis][selected_idx],
                        color="k",
                        linewidth=0.5,
                    )
                elif selected_axis == 1:
                    im.axes.axvline(
                        fine_idx[selected_axis][selected_idx],
                        color="k",
                        linewidth=0.5,
                    )
                plt.pause(0.1)
                fig.canvas.draw()
            # array of indices to keep
            keep_idx = np.concatenate(
                (
                    np.arange(selected_idx),
                    np.arange(
                        selected_idx + 1, coarse_axes[selected_axis].size
                    ),
                )
            )
            # renmove gridpoint
            coarse_axes[selected_axis] = coarse_axes[selected_axis][keep_idx]
            coarse_array = coarse_array.take(keep_idx, axis=selected_axis)
            # construct indices where new interpolated values belong
            # https://stackoverflow.com/a/42657219/5101335
            idx = [slice(None)] * fine_interped.ndim
            idx[selected_axis] = selected_interped_fine_idx
            # compute new interpolated
            fine_interped[idx] = selected_interped
            # compute residuals
            residuals = fine_interped - fine_array
            local_abs_rel_errs = np.abs(residuals / fine_array)
            max_local_abs_rel_err_idx = np.unravel_index(
                local_abs_rel_errs.argmax(), local_abs_rel_errs.shape
            )
            max_local_abs_rel_err = local_abs_rel_errs[
                max_local_abs_rel_err_idx
            ]
            max_local_abs_err = np.abs(residuals[max_local_abs_rel_err_idx])
            local_sqr_errs = (residuals * residuals) / fine_array.size
            global_err = local_sqr_errs.sum()
            fine_idx[selected_axis] = fine_idx[selected_axis][keep_idx]
            diagnostics = OrderedDict()
            diagnostics["removed axis".rjust(15)] = f"{selected_axis:15}"
            diagnostics["removed index".rjust(15)] = f"{selected_idx:15}"
            diagnostics[
                "max local abs. rel. err.".rjust(26)
            ] = f"{max_local_abs_rel_err:26.5E}"
            diagnostics[
                "corresponding abs. err.".rjust(25)
            ] = f"{max_local_abs_err:25.5E}"
            diagnostics["global err.".rjust(13)] = f"{global_err:13.5E}"
            diagnostics[
                "coarse axes shape".rjust(20)
            ] = f"{coarse_array.shape}".rjust(20)
            if len(removed) % 25 == 0:
                hlines = "".join(
                    ("-" * len(s.strip())).rjust(len(s))
                    for s in diagnostics.keys()
                )
                headings = "".join(diagnostics.keys())
                print("\n".join((hlines, headings, hlines)))
            print("".join(diagnostics.values()))
            removed.append((selected_axis, selected_idx))
        else:
            break
    with open("removed.pkl", "wb") as f:
        pickle.dump(removed, f)


def remove_removed(full_df, removed):
    fine_axes = [x.values for x in full_df.index.levels]
    fine_array = full_df.values.reshape(full_df.index.levshape)
    # remove removed indices
    coarse_axes = fine_axes.copy()
    coarse_array = fine_array
    for axis in range(fine_array.ndim):
        removed_indices = [index for a, index in removed if axis == a]
        coarse_array = np.delete(coarse_array, removed_indices, axis=axis)
        coarse_axes[axis] = np.delete(fine_axes[axis], removed_indices)
    interpolated = interpolate.RegularGridInterpolator(
        coarse_axes, coarse_array
    )(
        np.fromiter(
            itertools.chain(*itertools.product(*(fine_axes))), dtype=float
        ).reshape(fine_array.shape + (-1,))
    )


def split_into_monotonic_subsets(
    full_df,
    split_on="beta",
    rs=[5, 11],
    split_idx=217,
    value_at_max_cdf=1636.7475317348378,
    print_errors=False,
):
    # split ranges into subsets
    boundaries = [0, V.index.unique(split_on)[split_idx], np.inf]
    monotonic_subsets = []
    for i in range(len(rs)):
        # select a subset of the full dataset
        in_subset = (
            boundaries[i] <= full_df.columns.get_level_values(split_on)
        ) & (full_df.columns.get_level_values(split_on) < boundaries[i + 1])
        A_subset_full = full_df.loc[:, in_subset]
        # perform SVD on the subset
        U_subset, S_subset, Vt_subset = np.linalg.svd(
            A_subset_full, full_matrices=False
        )
        r = rs[i]
        U_subset_truncated = U_subset[:, :r]
        S_subset_truncated = S_subset[:r]
        Vt_subset_truncated = Vt_subset[:r, :]
        # reconstruct rank-r approximation
        A_subset_truncated = (
            U_subset_truncated * S_subset_truncated @ Vt_subset_truncated
        )
        # a monotonic entry is one one which is greater than or equal to all
        # entries that appeared previously
        monotonic_nondecreasing = (
            A_subset_truncated >= np.maximum.accumulate(A_subset_truncated)
        ).all(axis=1)
        # save monotonic subset
        monotonic_subset = pd.DataFrame(
            U_subset_truncated[monotonic_nondecreasing]
            * S_subset_truncated
            @ Vt_subset_truncated,
            index=U_df.index.unique("CDF")[monotonic_nondecreasing],
            columns=V_df.unstack().index[in_subset],
        )
        # check monoticity
        assert (monotonic_subset.diff().iloc[1:] > 0).all().all()
        monotonic_subsets.append(monotonic_subset)

    if print_errors:
        # linearly interpolate removed CDF values
        def interpolate_monotonic_subset(df):
            df = pd.concat(
                (
                    df,
                    pd.DataFrame(
                        np.nan,
                        index=full_df.index.unique("CDF").difference(df.index),
                        columns=df.columns,
                    ),
                )
            )
            # if the last CDF value is nonmonotonic the maximum value/CDF pair
            # will be necessary for interpolation
            if 1 not in df.index:
                df.loc[1] = value_at_max_cdf
                return df.sort_index().interpolate(method="index").iloc[:-1]
            else:
                return df.sort_index().interpolate(method="index")

        interpolated_monotonic_reconstruction = pd.concat(
            [interpolate_monotonic_subset(df) for df in monotonic_subsets],
            axis="columns",
        )
        # reconstruct the original matrix
        residuals = interpolated_monotonic_reconstruction - full_df
        l2_error = np.linalg.norm(residuals) / np.linalg.norm(full_df)
        linf_error = residuals.abs().max().max() / full_df.abs().max().max()
        print(
            f"rs: {rs}\n"
            f"split_idx: {split_idx}\n"
            f"l2_error: {l2_error}\n"
            f"linf_error: {linf_error}\n"
        )
    return monotonic_subsets


if __name__ == "__main__":
    # read full sab_data
    U_df = pd.read_hdf(
        "~/Developer/minimc/data/tsl/endfb8-fullorder-cdfrows/beta_endfb8_CDF_coeffs.hdf5"
    )
    S_df = pd.read_hdf(
        "~/Developer/minimc/data/tsl/endfb8-fullorder-cdfrows/beta_endfb8_S_coeffs.hdf5"
    )
    V_df = pd.read_hdf(
        "~/Developer/minimc/data/tsl/endfb8-fullorder-cdfrows/beta_endfb8_E_T_coeffs.hdf5"
    )
    # U_df = pd.read_hdf('~/Developer/minimc/data/tsl/endfb8-fullorder-cdfrows/alpha_endfb8_CDF_coeffs.hdf5')
    # S_df = pd.read_hdf('~/Developer/minimc/data/tsl/endfb8-fullorder-cdfrows/alpha_endfb8_S_coeffs.hdf5')
    # V_df = pd.read_hdf('~/Developer/minimc/data/tsl/endfb8-fullorder-cdfrows/alpha_endfb8_beta_T_coeffs.hdf5')
    # reconstruct full-rank matrix
    U = U_df["coefficient"].unstack()
    S = pd.DataFrame(
        np.diag(S_df.values.flatten()),
        index=S_df.index.get_level_values("order"),
        columns=S_df.index.get_level_values("order"),
    )
    V = V_df["coefficient"].unstack()
    A_full = U @ S @ V.T
    # split
    monotonic_subsets = split_into_monotonic_subsets(
        A_full, split_on="E", rs=[3, 5], split_idx=600, value_at_max_cdf=20
    )
    # res = split_subsets(A_full, split_on='beta', rs=[5,11], split_idx=217)
