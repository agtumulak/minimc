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
import matplotlib.pyplot as plt
import re
import tikzplotlib
from dask import dataframe as dd
from dask.diagnostics import ProgressBar
from functools import partial
from inspect import signature
from multiprocessing import Pool
from scipy import optimize, interpolate
from tqdm import tqdm

pbar = ProgressBar()
pbar.register()


# target nuclide mass ratio
A = 0.999167339

# Boltzmann constant in eV / K
k = 8.617333262145E-5

# bound cross section
sigma_free = 4.095600e1 / 2 # 2 hydrogens per H2O
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
        pattern = re.compile(r'\d([+-])')
        return float(pattern.sub(r'E\1', endf_float_string))

    with open(mf7_path, mode='r') as sab_file:
        # skip headers
        while True:
            if sab_file.readline().endswith('1 7  4\n'):
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
            for row in range(N_full_rows + (remainder != 0)):
                # Everything after column 66 is ignored
                line = sab_file.readline()
                doubles = [
                        to_float(line[start:start+11])
                        for start in range(0, 66, 11)
                        if not line[start:start+11].isspace()]
                for alpha, S in zip(doubles[::2], doubles[1::2]):
                    alphas.append(alpha)
                    betas.append(beta)
                    Ts.append(temp)
                    Ss.append(S)
            # The remaining temperatures are handled uniformly
            N_full_rows, remainder = divmod(N_alpha, 6)
            for _ in range(N_temp - 1):
                temp, beta = (
                        to_float(x) for x in sab_file.readline().split()[:2])
                # Subsequent betas use the first beta's alpha grid.
                unique_alphas = (a for a in alphas[:N_alpha])
                for row in range(N_full_rows + (remainder != 0)):
                    line = sab_file.readline()[:66]
                    for S in [
                            to_float(line[start:start+11])
                            for start in range(0, 66, 11)
                            if not line[start:start+11].isspace()]:
                        alphas.append(next(unique_alphas))
                        betas.append(beta)
                        Ts.append(temp)
                        Ss.append(S)
        df = pd.DataFrame.from_dict(
                {'alpha': alphas, 'beta': betas, 'T': Ts, 'S': Ss})
        # rescale alpha and beta due to LAT flag in LEAPR
        df['alpha'] = df['alpha'] * 293.6 / df['T']
        df['beta'] = df['beta'] * 293.6 / df['T']
        df = df.set_index(['beta', 'alpha', 'T'])
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
            numerators, denominators, out=numerators, where=denominators!=0)
    ratios = np.nan_to_num(ratios)
    return pd.Series(
            np.concatenate(([0], np.cumsum(ratios * (x[1:] - x[:-1])))),
            index=s.index)


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
    T_betas = np.array(sorted(sab_df.loc[:,:,T].index.unique('beta')))
    # valid beta values
    min_beta = - E / (k * T)
    max_beta = 20
    min_beta_index = np.searchsorted(T_betas, -min_beta) # as range end, will never include min_beta
    max_beta_index = np.searchsorted(T_betas, max_beta) # as range end, will never include max_beta
    betas = np.concatenate((
        [min_beta],
        -np.flip(T_betas[1:min_beta_index]),
        T_betas[:max_beta_index],
        [max_beta]))
    alpha_cdfs = {}
    alpha_pdfs = {}
    beta_pdf = pd.Series(0, index=betas) # pdf is zero at min_beta and max_beta
    for beta in betas[1:-1]:
        # energy-independent alpha distribution
        S_values = sab_df.xs((beta if beta > 0 else -beta, T), level=('beta', 'T'))['S']
        # energy-dependent alpha distribution, valid alpha values
        min_alpha = np.square(np.sqrt(E) - np.sqrt(E + beta * k * T)) / (A * k * T)
        max_alpha = np.square(np.sqrt(E) + np.sqrt(E + beta * k * T)) / (A * k * T)
        S_at_min_alpha, S_at_max_alpha = np.interp(
                [min_alpha, max_alpha], S_values.index, S_values)
        min_alpha_index = np.searchsorted(S_values.index, min_alpha, side='right') # will never include min_alpha
        max_alpha_index = np.searchsorted(S_values.index, max_alpha) # as range end, will never include max_alpha
        S_values = pd.concat((
                pd.Series({min_alpha: S_at_min_alpha}),
                S_values.iloc[min_alpha_index:max_alpha_index],
                pd.Series({max_alpha: S_at_max_alpha})))
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
    sab_s.index = sab_s.index.droplevel(['T', 'beta'])
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
    return c0 + c1 / T ** (1/2) + c2 / T + c3 / T ** (3/2) + c4 / T ** 2 + c5 / T ** (5/2)


def alpha_fitting_function(T, c0, c1, c2):
    return c0 + c1 / T + c2 / T ** 2


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
    label = label if label else 'direct'
    beta_pdf, alpha_pdfs, _, _, _ = process_E_T((sab_df, E, T))
    for beta, p_beta in beta_pdf.iloc[1:-1].iteritems():
        alpha_pdfs[beta] *= p_beta
    return (
            pd
            .concat(alpha_pdfs.values(), names=alpha_pdfs.keys(), axis='columns')
            .interpolate(method='index', axis='index', limit_area='inside')
            .fillna(0)
            .stack()
            .rename_axis(['alpha', 'beta'])
            .rename(label))


def get_pdf_pod(beta_T_path, beta_S_path, beta_E_CDF_path, alpha_T_path,
        alpha_S_path, alpha_beta_CDF_path, E, T, max_alpha, label=None):
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
    label = label if label else 'direct (pod reconstructed)'
    # Load data
    beta_T = pd.read_hdf(beta_T_path)
    beta_S = pd.read_hdf(beta_S_path)
    beta_E_CDF = pd.read_hdf(beta_E_CDF_path)
    alpha_T = pd.read_hdf(alpha_T_path)
    alpha_S = pd.read_hdf(alpha_S_path)
    alpha_beta_CDF = pd.read_hdf(alpha_beta_CDF_path)
    # evaluate beta at nearest E and nearest T
    Ts = beta_T.index.unique('T')
    nearest_T = Ts[np.argmin(np.abs(Ts - T))]
    Es = beta_E_CDF.index.unique('E')
    nearest_E = Es[np.argmin(np.abs(Es - E * 1e-6))]
    # Reconstruct beta CDF
    beta_cdf = pd.Series(
            (
                beta_T.loc[nearest_T].T.values
              @ np.diag(beta_S.values.flatten())
              @ beta_E_CDF.loc[nearest_E].unstack().T.values).flatten(),
            index=beta_E_CDF.index.unique('CDF'))
    # compute PDF
    beta_pdf = pd.Series(
            (beta_cdf.index[1:] - beta_cdf.index[:-1])
          / (beta_cdf.values[1:] - beta_cdf.values[:-1]),
            index=beta_cdf.values[:-1])
    # evaluate alpha at (possibly different) nearest T and nearest beta
    Ts = alpha_T.index.unique('T')
    nearest_T = Ts[np.argmin(np.abs(Ts - T))]
    # Reconstruct alpha CDF
    alpha_cdf = (
            pd.Series(
                (
                    alpha_T.loc[nearest_T].T.values
                  @ np.diag(alpha_S.values.flatten())
                  @ alpha_beta_CDF.unstack().T.values)
                .flatten(),
                index=alpha_beta_CDF.unstack().index)
            .unstack('beta'))
    # append alpha CDFs for negative beta values
    alpha_betas = alpha_cdf.columns
    # find largest beta in alpha_betas which is strictly less than E / (k * T)
    # we assume beta = 0 exists so result of searchsorted is >= 1
    min_beta = alpha_betas[np.searchsorted(alpha_betas, -beta_cdf.iloc[0]) - 1]
    neg_b_alpha_cdf = (
            alpha_cdf
            .loc[:, alpha_betas[1]:min_beta] # don't include beta = 0
            .rename(columns=lambda x: -x) # make beta labels negative
            .sort_index(axis='columns'))
    # find largest beta in alpha_betas which is strictly less than 20
    max_beta = alpha_betas[np.searchsorted(alpha_betas, 20) - 1]
    alpha_cdf = pd.concat(
            (
                neg_b_alpha_cdf,
                alpha_cdf.loc[:max_beta]),
            axis='columns')
    # add endpoints
    alpha_cdf.loc[0,:] = 0
    alpha_cdf.loc[1,:] = max_alpha
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
        s = pd.Series(s.index, index=pd.Index(s, name='alpha'))
        # insert values for min_alpha and max_apha
        min_alpha = np.square(np.sqrt(E) - np.sqrt(E + beta * k * T)) / (A * k * T)
        max_alpha = np.square(np.sqrt(E) + np.sqrt(E + beta * k * T)) / (A * k * T)
        s.loc[min_alpha] = np.nan # interpolate value later
        s.loc[max_alpha] = np.nan # interpolate value later

        s = s.sort_index().interpolate(method='index').loc[min_alpha:max_alpha]
        # rescale CDF to be 0 at min_alpha and 1 at max_alpha
        s = (s - s.min()) / (s.max() - s.min())
        alpha_pdf = pd.Series(
                (s.values[1:] - s.values[:-1]) / (s.index[1:] - s.index[:-1]),
                index=s.index[:-1])
        alpha_pdf.loc[min_alpha] = 0
        alpha_pdf.loc[max_alpha] = 0
        alpha_pdf = alpha_pdf.sort_index()
        # use common alpha grid
        new_alphas = set(alphas).difference(set(alpha_pdf.index))
        alpha_pdf = (
                pd.concat([alpha_pdf, pd.Series(np.nan, index=new_alphas)])
                .sort_index()
                .interpolate(method='index', limit_area='inside')
                .fillna(0)
                .loc[alphas])
        # interpolate value of beta pdf
        nonlocal beta_pdf
        beta_pdf_value = np.interp(beta, beta_pdf.index, beta_pdf, left=0, right=0)
        return alpha_pdf * beta_pdf_value
    bivariate_pdf = alpha_cdf.apply(get_joint_probability).stack()
    bivariate_pdf.index.names = ['alpha', 'beta']
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
    bin_areas = np.einsum('i,j->ij', *widths).reshape(-1)
    density = counts / bin_areas
    return pd.Series(
            density,
            index=pd.MultiIndex.from_product(
                bin_edges, names=['alpha', 'beta']),
            name='minimc')


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
    label = label if label else 'mcnp'
    cosine_bounds = [-1.,]
    energy_bounds = [0.,]
    counts = []
    with open(mctal_path) as f:
        line = ''
        # skip to cosine bins
        while not line.startswith('c'):
            line = f.readline()
        line = f.readline()
        # parse cosine boundaries until energy boundaries section is reached
        while not line.startswith('e'):
            cosine_bounds.extend(float(x) for x in line.split())
            line = f.readline()
        line = f.readline()
        # parse energy boundaries until time boundaries section is reached
        while not line.startswith('t'):
            energy_bounds.extend(float(x) * 1e6 for x in line.split())
            line = f.readline()
        # skip to values
        while not line.startswith('vals'):
            line = f.readline()
        line = f.readline()
        # parse values until the tfc section is reached
        while not line.startswith('tfc'):
            counts.extend(float(x) for x in line.split()[::2])
            line = f.readline()
    # compute densities in cosine, eV space
    cosine_bounds = np.array(cosine_bounds)
    energy_bounds = np.array(energy_bounds)
    counts = np.array(counts)
    cosine_widths = cosine_bounds[1:] - cosine_bounds[:-1]
    energy_widths = energy_bounds[1:] - energy_bounds[:-1]
    bin_areas = np.einsum('i,j->ij', cosine_widths, energy_widths).reshape(-1)
    density = counts / bin_areas
    s = pd.Series(
            density,
            index=pd.MultiIndex.from_product(
                [cosine_bounds[:-1], energy_bounds[:-1]], names=['mu', 'E']),
            name='mcnp')
    def to_alpha_beta(s):
        beta = s.name
        out_E = E + beta * k * T
        mus = s.index.get_level_values('mu')
        alphas = (E + out_E - 2 * mus * np.sqrt(E * out_E)) / (A * k * T)
        return pd.Series(s.values, index=alphas).rename_axis(index={'mu': 'alpha'})
    # multiply by jacobian
    s = s * 0.5 * A * (k * T) ** 2 / (np.sqrt(s.index.get_level_values('E') * E))
    s = (
            s
            .rename(index=lambda out_E: (out_E - E) / (k * T), level='E')
            .rename_axis(index={'E': 'beta'})
            .groupby('beta')
            .apply(to_alpha_beta))
    return (
            s[~s.index.duplicated()]
            .unstack('beta')
            .interpolate(method='index', axis='index', limit_area='inside')
            .fillna(0)
            .stack()
            .rename_axis(['alpha', 'beta'])
            .rename(label))


def get_pdf_minimc(minimc_path, E, T, label=None):
    """
    Creates bivariate PDF in alpha and beta from minimc output
    """
    label = label if label else 'minimc'
    with open(minimc_path) as f:
        line = ''
        # skip to cosine bins
        while not line.startswith('cosine'):
            line = f.readline()
        line = f.readline()
        cosine_bounds = [float(x) for x in f.readline().split(',')[:-1]]
        # skip to energy bins
        while not line.startswith('energy'):
            line = f.readline()
        line = f.readline()
        energy_bounds = [float(x) * 1e6 for x in f.readline().split(',')[:-1]]
        # skip to values
        while not line.startswith('mean'):
            line = f.readline()
        line = f.readline()
        counts = [float(x) for x in f.readline().split(',')[:-1]]
    # compute densities in cosine, eV space
    cosine_bounds = np.concatenate(([-np.inf], cosine_bounds, [np.inf]))
    energy_bounds = np.concatenate(([-np.inf], energy_bounds, [np.inf]))
    counts = np.array(counts)
    cosine_widths = cosine_bounds[1:] - cosine_bounds[:-1]
    energy_widths = energy_bounds[1:] - energy_bounds[:-1]
    bin_areas = np.einsum('i,j->ij', cosine_widths, energy_widths).reshape(-1)
    # compute densities and remove infinitely sized bins
    density = (
             (counts / bin_areas)
            .reshape(cosine_widths.size, energy_widths.size)
            [1:-1, 1:-1]
            .reshape(-1))
    # create pd.Series
    cosine_midpoints = (cosine_bounds[2:-1] + cosine_bounds[1:-2]) / 2
    energy_midpoints = (energy_bounds[2:-1] + energy_bounds[1:-2]) / 2
    s = pd.Series(
            density,
            index=pd.MultiIndex.from_product(
                [cosine_midpoints, energy_midpoints], names=['mu', 'E']),
            name='minimc')
    # convert to alpha, beta space; #TODO: this is identical to mcnp version
    def to_alpha_beta(s):
        beta = s.name
        out_E = E + beta * k * T
        mus = s.index.get_level_values('mu')
        alphas = (E + out_E - 2 * mus * np.sqrt(E * out_E)) / (A * k * T)
        return pd.Series(s.values, index=alphas).rename_axis(index={'mu': 'alpha'})
    # multiply by jacobian
    s = s * 0.5 * A * (k * T) ** 2 / (np.sqrt(s.index.get_level_values('E') * E))
    s = (
            s
            .rename(index=lambda out_E: (out_E - E) / (k * T), level='E')
            .rename_axis(index={'E': 'beta'})
            .groupby('beta')
            .apply(to_alpha_beta))
    return (
            s[~s.index.duplicated()]
            .unstack('beta')
            .interpolate(method='index', axis='index', limit_area='inside')
            .fillna(0)
            .stack()
            .rename_axis(['alpha', 'beta'])
            .rename(label))


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
            s.index.unique('beta'),
            s.index.unique('alpha'),
            np.log(s).unstack(),
            levels=100)
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$\alpha$')
    plt.title(r'$\log p_{\alpha, \beta} (\alpha, \beta)$')
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
    min_log_density = np.log(min(s[s>0].min() for s in series))
    min_beta = max(s.index.unique('beta').min() for s in series)
    max_beta = min(s.index.unique('beta').max() for s in series)
    min_alpha = max(s.index.unique('alpha').min() for s in series)
    max_alpha = min(s.index.unique('alpha').max() for s in series)
    f, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex='col', sharey='row', constrained_layout=True)
    for i, row in enumerate(axs):
        for j, ax in enumerate(row):
            k = ncols * i + j
            if k >= len(series):
                break
            s = series[k]
            cm = ax.contourf(
                s.index.unique('beta'),
                s.index.unique('alpha'),
                np.log(s).unstack(),
                levels=np.linspace(min_log_density, 0, 100))
            if (j == 0):
                ax.set_ylabel(r'$\alpha$')
            if (i == 1):
                ax.set_xlabel(r'$\beta$')
            ax.set_xlim(min_beta, max_beta)
            ax.set_ylim(min_alpha, max_alpha)
            ax.set_title(s.name)
    f.colorbar(cm, ax=axs, location='right')
    plt.show()


def marginalize(s, axis='beta'):
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
    return (
            s.groupby(axis)

            .apply(lambda s: np.trapz(s.values, s.index.droplevel(axis))))


def compare_univariate_pdf(title, *series, axis='beta'):
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
    plt.xlabel(fr'$\{axis}$')
    plt.ylabel(fr'$p_{{\{axis}}} (\{axis})$')
    plt.xlim(
            max(s.index.min() for s in series),
            min(s.index.max() for s in series))
    plt.legend()
    plt.title(title)
    plt.show()


def parallel_apply(df_grouped, func, *args):
    func_args = ((group_index, group) + args for group_index, group in df_grouped)
    with Pool(processes=8) as pool:
        return pd.concat(
                [x for x in tqdm(
                    pool.imap(func=func, iterable=func_args),
                    total=len(df_grouped))])


def fit_points(args):
    group_name, s = args
    s = pd.Series(
            coeffs := optimize.curve_fit(
                beta_fitting_function, s.index.unique(level='T'), s)[0],
            index=pd.Index(range(len(coeffs)), name='coefficient'))
    return pd.concat([s], keys=[group_name], names=['E', 'CDF'])


def beta_functional_expansion(sab_df, E_min=1e-5, E_max=4.0, n_Es=1000,
        n_cdfs=1000, order=3):
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
        Approximate number of CDF values to use
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
    df_Ts = np.array(sorted(sab_df.index.unique('T')))
    # equally spaced-lethargy intervals, do not include zero lethargy
    assert E_min > 0
    lethargies = np.linspace(0, np.log(E_max/E_min), num=n_Es+1)[:0:-1]
    Es = E_max * np.exp(-lethargies)
    with Pool(processes=8) as pool:
        # incident energy, temperature pairs
        E_T_values = np.array([(sab_df, E, T) for E in Es for T in df_Ts])
        results = (
                np.array(
                    [[x[2], x[4]] for x in tqdm(
                        pool.imap(func=process_E_T, iterable=E_T_values),
                        total=len(E_T_values))],
                    dtype=object).reshape(len(Es), len(df_Ts), -1))
        beta_cdfs = results[:,:,0]
        inelastic_xs = pd.DataFrame(results[:,:,1], index=Es, columns=df_Ts)
    F = np.linspace(0, 1, n_cdfs)
    beta_df = pd.DataFrame(
            np.nan,
            index=pd.Index(F, name='CDF'),
            columns=pd.MultiIndex.from_product((Es, df_Ts), names=('E', 'T')))
    # a good beta grid for the alpha functional expansion
    unique_betas = np.unique(np.abs(beta_df))
    N_betas = 1000
    selected_betas = unique_betas[
            np.round(np.linspace(0, unique_betas.size-1, N_betas)).astype(int)]
    # interpolate to fill in selected CDF values
    for E, x_E in zip(Es, beta_cdfs):
        for T, beta_cdf in zip(df_Ts, x_E):
            beta_df.loc[:, (E, T)] = np.interp(F, beta_cdf, beta_cdf.index)
    # perform proper orthogonal decomposition
    beta_df_pod_form = beta_df.stack('T').unstack('CDF')
    U, S, Vt = np.linalg.svd(beta_df_pod_form, full_matrices=False)
    order = S.size if order is None else order
    U_df = (pd.DataFrame(
        {'coefficient': U[:, :order].flatten()},
        index=pd.MultiIndex.from_product(
            [df_Ts, range(order)], names=['T', 'order'])))
    S_df = (pd.DataFrame(
        {'coefficient': S[:order]},
        index=pd.MultiIndex.from_product(
            [range(order)], names=['order'])))
    V_df = (pd.DataFrame(
        {'coefficient': Vt.T[:, :order].flatten()},
        index=pd.MultiIndex.from_product(
            [Es, F, range(order)], names=['E', 'CDF', 'order'])))
    # set energy units to MeV
    V_df.index = V_df.index.set_levels(
            V_df.index.unique(level='E') * 1e-6, level='E')
    # reconstruct
    beta_df_reconstructed = pd.DataFrame(
            U[:, :order] @ np.diag(S[:order]) @ Vt[:order, :],
            index=beta_df_pod_form.index,
            columns=beta_df_pod_form.columns)
    # check that CDFS are monotonic for certain T values
    print(f"RMSE: {np.sqrt(((beta_df_reconstructed - beta_df_pod_form)**2).mean().mean())}")
    monotonic_check_df = beta_df_reconstructed.stack('CDF').unstack('T')
    is_monotonic = monotonic_check_df.apply(lambda s: s.is_monotonic)
    print(
            f"Of {Es.size} incident energies and {df_Ts.size} target "
            f"temperatures, {is_monotonic.sum()} of {is_monotonic.size} "
            f"({is_monotonic.sum() / is_monotonic.size * 100}%) have "
            f"monotonic beta as a function of CDF")
    if not is_monotonic.all():
        print("The following CDFs are not monotonic:")
        print(monotonic_check_df.loc[:, ~is_monotonic])
    return U_df, S_df, V_df


def fit_alpha(args):
    """
    Fits alpha as a function of temperature for a given beta and CDF value.

    Parameters
    ----------
    args : tuple of (group key, group)
        The first argument is the group key (a (beta, CDF) pair).
        The second element is a pd.Series containing a MultiIndex. The
        MultiIndex levels are `beta`, `CDF`, and temperature `T`. There is only
        a single value of `beta` and `CDF` while multiple `T` values must be
        present. The values are values of alpha.
    """
    [beta, F], s= args
    try:
        coeffs = optimize.curve_fit(
                alpha_fitting_function, s.index.get_level_values('T'), s)[0]
    except TypeError:
        # This happens when there are N coefficients to fit, M data points, and
        # N > M. In this case, we only fit the first M coefficients and set the
        # other N - M coefficients to zero.
        n_coeffs = len(signature(alpha_fitting_function).parameters) - 1
        kwargs = {f'c{i}': 0 for i in range(len(s), n_coeffs)}
        coeffs = optimize.curve_fit(
                partial(alpha_fitting_function, **kwargs),
                s.index.get_level_values('T'), s)[0]
        coeffs = np.concatenate([coeffs, np.zeros(n_coeffs - len(s))])
    return pd.Series(
            coeffs,
            index=pd.MultiIndex.from_product(
                [[beta], [F], range(len(coeffs))],
                names=['beta', 'CDF', 'coefficient']))


def alpha_functional_expansion(sab_df, selected_betas, n_cdfs=1000, order=3):
    """
    Computes the energy-independent conditional CDF in alpha given beta at
    various beta values and temperatures, then performs a functional expansion
    in temperature at various beta values and CDF values.

    Parameters
    ----------
    sab_df : pd.DataFrame
        S(a,b,T) DataFrame
    n_betas : np.ndarray
        Beta values to use
    n_cdfs : int, optional
        Approximate number of CDF values to use
    order : int, optional
        Expansion order for proper orthogonal decomposition. Setting to None
        will return the full expansion.
    """
    df_Ts = np.array(sorted(sab_df.index.unique('T')))
    selected_betas = set(selected_betas)
    def common_beta_grid(group):
        """
        Modifies beta values at this temperature to conform with common beta
        grid.
        """
        nonlocal sab_df
        T = group.index.unique('T')[0]
        group.index = group.index.droplevel('T')
        # add betas from selected_betas which are not already in group
        new_betas = selected_betas.difference(group.index.unique('beta'))
        new_df = pd.DataFrame(
                np.nan,
                columns=group.index.unique('alpha'),
                index=pd.Index(new_betas, name='beta'))
        combined_df = pd.concat([group.unstack('alpha')['S'], new_df]).sort_index()
        # S(a,b) above maximum beta at this temperature is zero
        combined_df.loc[combined_df.index > group.index.unique('beta').max()] = 0
        # interpolate linearly in beta and linearly in ln(S) (interpolation
        # scheme 4 in ENDF)
        return (
                np.exp(
                    np.log(combined_df)
                    .interpolate(method='index', axis='index', limit_area='inside'))
                .loc[selected_betas]
                .stack())
    common_beta_sab_df = sab_df.groupby('T').apply(common_beta_grid)
    # alpha grids do not have to match across T-beta pairs, but they do have to
    # have the same minimum and maximum values
    print("computing alpha CDFs...")
    largest_alpha = sab_df.index.unique('alpha').max()
    print(f"largest alpha: {largest_alpha}")
    alpha_cdfs = common_beta_sab_df.groupby(['T', 'beta']).apply(process_b_T, largest_alpha)
    # take the union of all CDF values that appear across all incident energies
    all_cdfs = sorted(set(alpha_cdfs))
    # choose number of CDF points we want to use
    F = all_cdfs[::len(all_cdfs) // n_cdfs]
    # don't include 0
    F = F if F[0] != 0 else F[1:]
    # last CDF must always be 1.
    print(f"using {len(F)} CDF values...")
    # interpolate alpha values at selected CDF values
    alpha_df = alpha_cdfs.groupby(['T', 'beta']).apply(
            lambda s:
            pd.Series(
                np.interp(F, s, s.index.get_level_values('alpha')),
                index=pd.Index(F, name='CDF')))
    # construct matrix to decompose with POD
    alpha_df_pod_form = alpha_df.unstack(['beta', 'CDF'])
    # impute NaN entries: https://stackoverflow.com/a/35611142/5101335
    nan_entries = alpha_df_pod_form.isna()
    alpha_df_pod_form = alpha_df_pod_form.fillna(value=alpha_df_pod_form.mean())
    # repeat multiple times until convergence
    converged = False
    prev_trace_norm = np.nan
    while not converged:
        U, S, Vt = np.linalg.svd(alpha_df_pod_form, full_matrices=False)
        order = S.size if order is None else order
        A = U[:, :order] @ np.diag(S[:order]) @ Vt[:order,:]
        alpha_df_pod_form = alpha_df_pod_form.where(~nan_entries, A)
        trace_norm = S.sum()
        abs_rel_diff = np.abs((trace_norm - prev_trace_norm) / prev_trace_norm)
        print (f"trace norm abs rel diff: {abs_rel_diff}", end=" " *  10 + '\r')
        if abs_rel_diff <  1e-8:
            break
        else:
            prev_trace_norm = trace_norm
    U_df = (pd.DataFrame(
        {'coefficient': U[:, :order].flatten()},
        index=pd.MultiIndex.from_product(
            [df_Ts, range(order)], names=['T', 'order'])))
    S_df = (pd.DataFrame(
        {'coefficient': S[:order]},
        index=pd.MultiIndex.from_product(
            [range(order)], names=['order'])))
    V_df = (pd.DataFrame(
        {'coefficient': Vt.T[:, :order].flatten()},
        index=pd.MultiIndex.from_product([
            alpha_df.index.unique('beta'), F, range(order)],
            names=['beta', 'CDF', 'order'])))
    # reconstruct
    alpha_df_reconstructed = pd.DataFrame(
            U[:, :order] @ np.diag(S[:order]) @ Vt[:order, :],
            index=alpha_df_pod_form.index,
            columns=alpha_df_pod_form.columns)
    print(f"RMSE: {np.sqrt(((alpha_df_reconstructed - alpha_df_pod_form)**2).mean().mean())}")
    # check that CDFS are monotonic for certain T values
    monotonic_check_df = alpha_df_reconstructed.stack('CDF').unstack('T')
    monotonic_check_df = monotonic_check_df.loc[:, ~monotonic_check_df.isna().any()]
    is_monotonic = monotonic_check_df.apply(lambda s: s.is_monotonic)
    print(
            f"Of {alpha_df.index.unique('beta').size} betas and {df_Ts.size} "
            f"target temperatures, there are {is_monotonic.size} non-NaN "
            f"columns. {is_monotonic.sum()} of {is_monotonic.size} "
            f"({is_monotonic.sum() / is_monotonic.size * 100}%) have "
            f"monotonic alpha as a function of CDF")
    if not is_monotonic.all():
        print("The following CDFs are not monotonic:")
        print(monotonic_check_df.loc[:, ~is_monotonic])
    return U_df, S_df, V_df


def uniform_coarsen(U, S, V):
    # helper function
    def compute_rel_frobenius_norm(row, full_df, interpolation_points):
        idx = pd.IndexSlice
        fine_Ts = full_df.index
        fine_betas = full_df.columns.unique('beta')
        fine_CDFs = full_df.columns.unique('CDF')
        coarse_Ts = fine_Ts[np.round(np.linspace(0, fine_Ts.size - 1, row['T gridsize'])).astype(int)]
        coarse_betas = fine_betas[np.round(np.linspace(0, fine_betas.size - 1, row['beta gridsize'])).astype(int)]
        coarse_CDFs = fine_CDFs[np.round(np.linspace(0, fine_CDFs.size - 1, row['CDF gridsize'])).astype(int)]
        coarse_interpolator = interpolate.RegularGridInterpolator(
                (coarse_Ts, coarse_betas, coarse_CDFs),
                (
                    full_df
                    .loc[coarse_Ts, idx[coarse_betas, coarse_CDFs]]
                    .values
                    .reshape(coarse_Ts.size, coarse_betas.size, coarse_CDFs.size)))
        interpolated_df = pd.DataFrame(
                coarse_interpolator(interpolation_points),
                index=full_df.index,
                columns=full_df.columns)
        norm = np.linalg.norm(full_df - interpolated_df) / np.linalg.norm(full_df)
        return norm
    # parse hdf5 data
    U = U.unstack()
    S = pd.DataFrame(np.diag(S.values.flatten()), index=U.columns, columns=U.columns)
    V = V.unstack()
    full_df = U @ S @ V.T
    # get full grids
    fine_Ts = U.index
    fine_betas = V.index.unique('beta')
    fine_CDFs = V.index.unique('CDF')
    # points which will be interpolated
    interpolation_points = (
            np.fromiter(
                (
                    val
                    for T in fine_Ts
                    for beta in fine_betas
                    for CDF in fine_CDFs
                    for val in (T, beta, CDF)),
                float,
                count=fine_Ts.size * fine_betas.size * fine_CDFs.size * 3)
            .reshape(fine_Ts.size, fine_betas.size * fine_CDFs.size, -1))
    # get coarsened gridsizes
    T_gridsizes = range(2, fine_Ts.size + 1)
    E_gridsizes = range(100, fine_betas.size + 1, 100)
    CDF_gridsizes = range(100, fine_CDFs.size + 1, 100)
    # initialize Series of relative errors
    index = pd.MultiIndex.from_product(
            (T_gridsizes, E_gridsizes, CDF_gridsizes),
            names=['T gridsize', 'beta gridsize', 'CDF gridsize'])
    rel_frobenius_norms = (
            pd.Series(0, index=index, name='relative error')
            .reset_index())
    ddf = dd.from_pandas(rel_frobenius_norms, npartitions=8)
    ddf['relative error'] = ddf.apply(
            compute_rel_frobenius_norm,
            axis='columns',
            args=(full_df, interpolation_points), meta='double')
    rel_frobenius_norms = ddf.compute().set_index(['T gridsize', 'beta gridsize', 'CDF gridsize'])
    return rel_frobenius_norms


if __name__ == '__main__':
    # Broomstick problem
    E = 0.56 # energy in eV
    T = 450.0 # temperature in K
    pdf_mcnp = get_pdf_mcnp(f'/Users/atumulak/Developer/mcnp6-runs/broomstick/endf80-{E}eV-{T:.1f}K.mctal', E, T, label='MCNP (tabular)')
    pdf_minimc = get_pdf_minimc(f'/Users/atumulak/Developer/minimc/broomstick_{E}eV_{T:.1f}K.minimc', E, T, label='MiniMC (POD)')

    plt.clf()
    plt.ticklabel_format(style='plain')
    beta_pdf_mcnp = marginalize(pdf_mcnp, axis='beta')
    beta_pdf_minimc = marginalize(pdf_minimc, axis='beta')
    beta_pdf_mcnp.plot(linestyle='-', color='k')
    beta_pdf_minimc.plot(linestyle=':', color='k')
    # plt.tight_layout()
    plt.ylabel(r'$p_{\beta}(\beta)$')
    plt.xlabel(r'$\beta$')
    plt.xlim(max(beta_pdf_mcnp.index.min(), beta_pdf_minimc.index.min()), 10)
    plt.ylim(0, 0.25)
    plt.legend(frameon=False, loc='upper left')
    plt.grid(True)
    tikzplotlib.clean_figure()
    tikzplotlib.save('beta_comparison.tex')

    plt.clf()
    alpha_pdf_mcnp = marginalize(pdf_mcnp, axis='alpha')
    alpha_pdf_minimc = marginalize(pdf_minimc, axis='alpha')
    alpha_pdf_mcnp.plot(linestyle='-', color='k')
    alpha_pdf_minimc.plot(linestyle=':', color='k')
    # plt.tight_layout()
    plt.ylabel(r'$p_{\alpha}(\alpha)$')
    plt.xlabel(r'$\alpha$')
    plt.xlim(0, 50)
    plt.ylim(0, 0.07)
    plt.legend(frameon=False)
    plt.grid(True)
    tikzplotlib.clean_figure()
    tikzplotlib.save('alpha_comparison.tex')


    # Segmented Problem

    # MCNP
    plt.clf()
    mcnp_x = np.array([0.0, 1.3888E-11, 1.4747E-11, 1.5659E-11, 1.6627E-11, 1.7655E-11, 1.8747E-11, 1.9906E-11, 2.1137E-11, 2.2444E-11, 2.3832E-11, 2.5305E-11, 2.6870E-11, 2.8532E-11, 3.0296E-11, 3.2170E-11, 3.4159E-11, 3.6271E-11, 3.8514E-11, 4.0896E-11, 4.3424E-11, 4.6110E-11, 4.8961E-11, 5.1988E-11, 5.5203E-11, 5.8617E-11, 6.2241E-11, 6.6090E-11, 7.0177E-11, 7.4517E-11, 7.9124E-11, 8.4017E-11, 8.9212E-11, 9.4729E-11, 1.0059E-10, 1.0681E-10, 1.1341E-10, 1.2042E-10, 1.2787E-10, 1.3578E-10, 1.4417E-10, 1.5309E-10, 1.6256E-10, 1.7261E-10, 1.8328E-10, 1.9461E-10, 2.0665E-10, 2.1943E-10, 2.3300E-10, 2.4740E-10, 2.6270E-10, 2.7895E-10, 2.9620E-10, 3.1451E-10, 3.3396E-10, 3.5461E-10, 3.7654E-10, 3.9982E-10, 4.2455E-10, 4.5080E-10, 4.7867E-10, 5.0827E-10, 5.3970E-10, 5.7308E-10, 6.0851E-10, 6.4614E-10, 6.8610E-10, 7.2852E-10, 7.7357E-10, 8.2141E-10, 8.7220E-10, 9.2614E-10, 9.8341E-10, 1.0442E-09, 1.1088E-09, 1.1774E-09, 1.2502E-09, 1.3275E-09, 1.4095E-09, 1.4967E-09, 1.5893E-09, 1.6875E-09, 1.7919E-09, 1.9027E-09, 2.0203E-09, 2.1453E-09, 2.2779E-09, 2.4188E-09, 2.5684E-09, 2.7272E-09, 2.8958E-09, 3.0749E-09, 3.2650E-09, 3.4669E-09, 3.6813E-09, 3.9089E-09, 4.1507E-09, 4.4073E-09, 4.6798E-09, 4.9692E-09, 5.2765E-09, 5.6028E-09, 5.9493E-09, 6.3171E-09, 6.7078E-09, 7.1225E-09, 7.5630E-09, 8.0307E-09, 8.5272E-09, 9.0545E-09, 9.6144E-09, 1.0209E-08, 1.0840E-08, 1.1511E-08, 1.2222E-08, 1.2978E-08, 1.3781E-08, 1.4633E-08, 1.5538E-08, 1.6498E-08, 1.7519E-08, 1.8602E-08, 1.9752E-08, 2.0974E-08, 2.2271E-08, 2.3648E-08, 2.5110E-08, 2.6663E-08, 2.8311E-08, 3.0062E-08, 3.1921E-08, 3.3895E-08, 3.5991E-08, 3.8216E-08, 4.0580E-08, 4.3089E-08, 4.5753E-08, 4.8583E-08, 5.1587E-08, 5.4777E-08, 5.8164E-08, 6.1761E-08, 6.5580E-08, 6.9635E-08, 7.3941E-08, 7.8513E-08, 8.3368E-08, 8.8523E-08, 9.3997E-08, 9.9810E-08, 1.0598E-07, 1.1254E-07, 1.1949E-07, 1.2688E-07, 1.3473E-07, 1.4306E-07, 1.5191E-07, 1.6130E-07, 1.7127E-07, 1.8187E-07, 1.9311E-07, 2.0505E-07, 2.1773E-07, 2.3120E-07, 2.4549E-07, 2.6067E-07, 2.7679E-07, 2.9391E-07, 3.1208E-07, 3.3138E-07, 3.5187E-07, 3.7363E-07, 3.9673E-07, 4.2127E-07, 4.4732E-07, 4.7498E-07, 5.0435E-07, 5.3553E-07, 5.6865E-07, 6.0381E-07, 6.4115E-07, 6.8080E-07, 7.2290E-07, 7.6760E-07, 8.1506E-07, 8.6546E-07, 9.1898E-07, 9.7581E-07, 1.0361E-06, 1.1002E-06, 1.1683E-06, 1.2405E-06, 1.3172E-06, 1.3987E-06, 1.4851E-06, 1.5770E-06, 1.6745E-06, 1.7780E-06, 1.8880E-06, 2.0047E-06, 2.1287E-06, 2.2603E-06])
    mcnp_y = np.array([0.00000E+00, 0.00000E+00, 0.00000E+00, 8.47025E-08, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 8.70354E-08, 6.19715E-08, 0.00000E+00, 0.00000E+00, 0.00000E+00, 9.16434E-08, 7.96543E-08, 0.00000E+00, 0.00000E+00, 7.24145E-08, 0.00000E+00, 4.21454E-08, 0.00000E+00, 1.18403E-07, 0.00000E+00, 6.70947E-08, 2.40811E-07, 0.00000E+00, 0.00000E+00, 1.39505E-07, 1.42467E-07, 1.65032E-07, 0.00000E+00, 5.71574E-08, 0.00000E+00, 1.53026E-07, 1.97758E-07, 2.42679E-07, 3.59412E-07, 4.34465E-07, 4.13759E-07, 3.32021E-07, 7.08945E-07, 7.88901E-07, 8.10570E-07, 8.38619E-07, 1.07238E-06, 7.97081E-07, 1.46928E-06, 2.08201E-06, 1.37004E-06, 1.26826E-06, 1.72541E-06, 2.04854E-06, 2.02357E-06, 2.00285E-06, 3.10952E-06, 3.54243E-06, 3.42757E-06, 5.25450E-06, 5.07330E-06, 5.11579E-06, 6.59894E-06, 8.46631E-06, 7.86158E-06, 9.81307E-06, 1.11746E-05, 1.26110E-05, 1.35622E-05, 1.42971E-05, 1.61058E-05, 2.06618E-05, 2.36974E-05, 2.54540E-05, 2.96002E-05, 3.60012E-05, 3.86273E-05, 4.49236E-05, 5.02828E-05, 5.45810E-05, 6.31841E-05, 7.45957E-05, 8.29121E-05, 9.26663E-05, 1.07870E-04, 1.22184E-04, 1.39470E-04, 1.55324E-04, 1.75970E-04, 1.94682E-04, 2.18780E-04, 2.55589E-04, 2.73570E-04, 3.15368E-04, 3.61919E-04, 3.99988E-04, 4.56209E-04, 5.12500E-04, 5.81657E-04, 6.62079E-04, 7.40609E-04, 8.26175E-04, 9.29266E-04, 1.04715E-03, 1.14894E-03, 1.29892E-03, 1.42747E-03, 1.61018E-03, 1.75669E-03, 1.98022E-03, 2.19032E-03, 2.42768E-03, 2.67316E-03, 2.90153E-03, 3.20034E-03, 3.47882E-03, 3.82311E-03, 4.14036E-03, 4.44038E-03, 4.84306E-03, 5.14048E-03, 5.52561E-03, 5.85232E-03, 6.26373E-03, 6.58475E-03, 6.88797E-03, 7.20599E-03, 7.46233E-03, 7.68745E-03, 7.92206E-03, 8.03783E-03, 8.07696E-03, 8.10170E-03, 8.07846E-03, 7.97689E-03, 7.79583E-03, 7.53001E-03, 7.19942E-03, 6.86645E-03, 6.47997E-03, 5.96922E-03, 5.52570E-03, 5.05106E-03, 4.50534E-03, 3.99768E-03, 3.56550E-03, 3.13483E-03, 2.68274E-03, 2.34219E-03, 2.01003E-03, 1.76915E-03, 1.51363E-03, 1.37651E-03, 1.21478E-03, 1.10748E-03, 1.04479E-03, 1.01496E-03, 1.00566E-03, 9.93991E-04, 2.80501E-03, 6.14868E-04, 2.75682E-04, 1.23904E-04, 5.47378E-05, 2.18376E-05, 8.40966E-06, 3.28966E-06, 1.36932E-06, 9.91850E-07, 2.98193E-07, 9.91244E-08, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00])
    mcnp_midpoints = (mcnp_x[1:] + mcnp_x[:-1]) / 2
    mcnp_segments_density = mcnp_y / (mcnp_x[1:] - mcnp_x[:-1])
    plt.plot(
            mcnp_midpoints * 1e6, mcnp_segments_density / 1e6, label='MCNP (tabular)',
            color='k', linestyle='-')

    # MINIMC
    pod_x = np.array([
        -np.inf,
        1.388794e-11, 1.474673e-11, 1.565861e-11, 1.662689e-11, 1.765504e-11, 1.874676e-11, 1.990600e-11, 2.113692e-11, 2.244395e-11, 2.383181e-11, 2.530548e-11, 2.687029e-11, 2.853185e-11, 3.029616e-11, 3.216957e-11, 3.415883e-11, 3.627109e-11, 3.851397e-11, 4.089554e-11, 4.342438e-11, 4.610960e-11, 4.896086e-11, 5.198843e-11, 5.520321e-11, 5.861679e-11, 6.224145e-11, 6.609024e-11, 7.017703e-11, 7.451654e-11, 7.912438e-11, 8.401716e-11, 8.921249e-11, 9.472909e-11, 1.005868e-10, 1.068067e-10, 1.134113e-10, 1.204243e-10, 1.278709e-10, 1.357780e-10, 1.441740e-10, 1.530893e-10, 1.625558e-10, 1.726077e-10, 1.832811e-10, 1.946146e-10, 2.066489e-10, 2.194273e-10, 2.329960e-10, 2.474036e-10, 2.627022e-10, 2.789468e-10, 2.961959e-10, 3.145116e-10, 3.339600e-10, 3.546109e-10, 3.765388e-10, 3.998227e-10, 4.245463e-10, 4.507988e-10, 4.786746e-10, 5.082742e-10, 5.397041e-10, 5.730776e-10, 6.085147e-10, 6.461432e-10, 6.860984e-10, 7.285244e-10, 7.735738e-10, 8.214090e-10, 8.722020e-10, 9.261360e-10, 9.834051e-10, 1.044215e-09, 1.108786e-09, 1.177350e-09, 1.250153e-09, 1.327458e-09, 1.409543e-09, 1.496705e-09, 1.589256e-09, 1.687530e-09, 1.791881e-09, 1.902685e-09, 2.020340e-09, 2.145271e-09, 2.277927e-09, 2.418786e-09, 2.568356e-09, 2.727174e-09, 2.895813e-09, 3.074880e-09, 3.265020e-09, 3.466917e-09, 3.681300e-09, 3.908938e-09, 4.150654e-09, 4.407316e-09, 4.679849e-09, 4.969235e-09, 5.276515e-09, 5.602796e-09, 5.949254e-09, 6.317135e-09, 6.707765e-09, 7.122550e-09, 7.562984e-09, 8.030653e-09, 8.527241e-09, 9.054536e-09, 9.614437e-09, 1.020896e-08, 1.084025e-08, 1.151057e-08, 1.222234e-08, 1.297813e-08, 1.378066e-08, 1.463280e-08, 1.553765e-08, 1.649844e-08, 1.751865e-08, 1.860194e-08, 1.975222e-08, 2.097363e-08, 2.227056e-08, 2.364770e-08, 2.510999e-08, 2.666271e-08, 2.831144e-08, 3.006212e-08, 3.192106e-08, 3.389494e-08, 3.599089e-08, 3.821644e-08, 4.057961e-08, 4.308892e-08, 4.575339e-08, 4.858262e-08, 5.158680e-08, 5.477675e-08, 5.816395e-08, 6.176061e-08, 6.557968e-08, 6.963490e-08, 7.394088e-08, 7.851313e-08, 8.336811e-08, 8.852330e-08, 9.399728e-08, 9.980975e-08, 1.059816e-07, 1.125352e-07, 1.194940e-07, 1.268831e-07, 1.347291e-07, 1.430602e-07, 1.519066e-07, 1.613000e-07, 1.712742e-07, 1.818652e-07, 1.931111e-07, 2.050525e-07, 2.177322e-07, 2.311960e-07, 2.454924e-07, 2.606728e-07, 2.767919e-07, 2.939077e-07, 3.120820e-07, 3.313800e-07, 3.518714e-07, 3.736299e-07, 3.967339e-07, 4.212666e-07, 4.473162e-07, 4.749767e-07, 5.043477e-07, 5.355348e-07, 5.686504e-07, 6.038138e-07, 6.411515e-07, 6.807981e-07, 7.228963e-07, 7.675977e-07, 8.150633e-07, 8.654640e-07, 9.189814e-07, 9.758080e-07, 1.036149e-06, 1.100220e-06, 1.168254e-06, 1.240495e-06, 1.317203e-06, 1.398654e-06, 1.485142e-06, 1.576978e-06, 1.674493e-06, 1.778038e-06, 1.887986e-06, 2.004732e-06, 2.128698e-06, 2.260329e-06,
        np.inf])
    pod_midpoints = (pod_x[2:-1] + pod_x[1:-2]) / 2

    # Piecewise constant temperature
    pod_segments_y = np.array([0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e-07, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e-07, 1.000000e-07, 0.000000e+00, 1.000000e-07, 4.000000e-07, 0.000000e+00, 1.000000e-07, 0.000000e+00, 0.000000e+00, 1.000000e-07, 2.000000e-07, 1.000000e-07, 1.000000e-07, 5.000000e-07, 2.000000e-07, 4.000000e-07, 4.000000e-07, 5.000000e-07, 5.000000e-07, 8.000000e-07, 6.000000e-07, 1.200000e-06, 7.000000e-07, 8.000000e-07, 1.400000e-06, 2.000000e-06, 1.600000e-06, 3.600000e-06, 2.100000e-06, 2.700000e-06, 3.600000e-06, 4.000000e-06, 4.700000e-06, 5.500000e-06, 6.500000e-06, 5.500000e-06, 8.000000e-06, 1.110000e-05, 1.270000e-05, 1.160000e-05, 1.350000e-05, 1.670000e-05, 1.760000e-05, 2.080000e-05, 2.350000e-05, 2.790000e-05, 3.050000e-05, 3.790000e-05, 4.040000e-05, 4.570000e-05, 5.330000e-05, 6.190000e-05, 6.870000e-05, 8.110000e-05, 8.650000e-05, 1.056000e-04, 1.153000e-04, 1.291000e-04, 1.489000e-04, 1.769000e-04, 1.929000e-04, 2.158000e-04, 2.406000e-04, 2.807000e-04, 3.118000e-04, 3.516000e-04, 3.917000e-04, 4.542000e-04, 5.108000e-04, 5.815000e-04, 6.491000e-04, 7.273000e-04, 8.277000e-04, 9.230000e-04, 1.038400e-03, 1.168800e-03, 1.311900e-03, 1.434400e-03, 1.619600e-03, 1.801400e-03, 2.012800e-03, 2.222600e-03, 2.469000e-03, 2.704100e-03, 2.929900e-03, 3.215500e-03, 3.527300e-03, 3.836800e-03, 4.133600e-03, 4.494900e-03, 4.857800e-03, 5.211200e-03, 5.558700e-03, 5.922800e-03, 6.348000e-03, 6.611400e-03, 6.971200e-03, 7.234000e-03, 7.481500e-03, 7.788100e-03, 7.958200e-03, 8.075300e-03, 8.181400e-03, 8.227500e-03, 8.190300e-03, 8.094600e-03, 7.858100e-03, 7.650200e-03, 7.333100e-03, 6.977800e-03, 6.527800e-03, 6.092900e-03, 5.583100e-03, 5.104100e-03, 4.581900e-03, 4.094900e-03, 3.632300e-03, 3.170700e-03, 2.737700e-03, 2.378900e-03, 2.062700e-03, 1.777600e-03, 1.556700e-03, 1.355200e-03, 1.230900e-03, 1.146500e-03, 1.074700e-03, 1.071000e-03, 1.072200e-03, 1.063700e-03, 2.912700e-03, 7.005000e-04, 3.262000e-04, 1.793000e-04, 1.080000e-04, 7.760000e-05, 5.850000e-05, 5.540000e-05, 4.530000e-05, 4.400000e-05, 3.320000e-05, 2.040000e-05, 1.190000e-05, 5.400000e-06, 3.800000e-06, 1.400000e-06, 1.000000e-06, 4.000000e-07, 2.000000e-07, 2.000000e-07, 1.000000e-07, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])
    pod_segments_density = pod_segments_y[1:-1] / (pod_x[2:-1] - pod_x[1:-2])
    plt.plot(
            pod_midpoints * 1e6, pod_segments_density / 1e6, label=r'MiniMC(POD, piecewise constant $T$)',
            color='k', linestyle=':')

    # Continuous temperature
    pod_continuous_y = np.array([0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e-07, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e-07, 0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e-07, 0.000000e+00, 1.000000e-07, 2.000000e-07, 0.000000e+00, 1.000000e-07, 1.000000e-07, 2.000000e-07, 2.000000e-07, 6.000000e-07, 3.000000e-07, 2.000000e-07, 6.000000e-07, 4.000000e-07, 1.100000e-06, 1.200000e-06, 8.000000e-07, 1.100000e-06, 1.800000e-06, 1.900000e-06, 2.800000e-06, 2.900000e-06, 2.700000e-06, 3.800000e-06, 4.100000e-06, 3.800000e-06, 5.700000e-06, 5.200000e-06, 5.700000e-06, 9.200000e-06, 8.900000e-06, 1.170000e-05, 1.380000e-05, 1.160000e-05, 1.550000e-05, 1.650000e-05, 2.300000e-05, 2.180000e-05, 2.340000e-05, 3.140000e-05, 3.350000e-05, 3.820000e-05, 4.250000e-05, 4.950000e-05, 5.320000e-05, 6.210000e-05, 6.900000e-05, 8.100000e-05, 9.150000e-05, 1.100000e-04, 1.207000e-04, 1.415000e-04, 1.614000e-04, 1.748000e-04, 1.994000e-04, 2.178000e-04, 2.507000e-04, 2.846000e-04, 3.281000e-04, 3.645000e-04, 4.061000e-04, 4.634000e-04, 5.353000e-04, 5.803000e-04, 6.718000e-04, 7.540000e-04, 8.341000e-04, 9.396000e-04, 1.059500e-03, 1.189400e-03, 1.333200e-03, 1.472600e-03, 1.658500e-03, 1.785000e-03, 2.026800e-03, 2.234100e-03, 2.479200e-03, 2.774400e-03, 3.001500e-03, 3.277100e-03, 3.581400e-03, 3.902000e-03, 4.251800e-03, 4.570600e-03, 4.958100e-03, 5.276900e-03, 5.709700e-03, 5.954100e-03, 6.341600e-03, 6.656100e-03, 6.967900e-03, 7.277500e-03, 7.543400e-03, 7.794500e-03, 7.951900e-03, 8.086900e-03, 8.148900e-03, 8.217200e-03, 8.150000e-03, 7.964300e-03, 7.807400e-03, 7.556800e-03, 7.300100e-03, 6.897300e-03, 6.427400e-03, 5.975200e-03, 5.465400e-03, 4.995700e-03, 4.501000e-03, 4.028800e-03, 3.523000e-03, 3.101200e-03, 2.684100e-03, 2.304300e-03, 1.997200e-03, 1.712900e-03, 1.522500e-03, 1.316900e-03, 1.193800e-03, 1.114000e-03, 1.100800e-03, 1.057800e-03, 1.080700e-03, 1.050000e-03, 2.885400e-03, 7.053000e-04, 3.313000e-04, 1.850000e-04, 1.047000e-04, 7.450000e-05, 6.060000e-05, 5.420000e-05, 4.730000e-05, 4.300000e-05, 3.260000e-05, 1.810000e-05, 8.600000e-06, 4.800000e-06, 2.600000e-06, 1.600000e-06, 7.000000e-07, 4.000000e-07, 1.000000e-07, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])
    pod_continuous_density = pod_continuous_y[1:-1] / (pod_x[2:-1] - pod_x[1:-2])
    plt.plot(
            pod_midpoints * 1e6, pod_continuous_density / 1e6, label=r'MiniMC(POD, linear $T$)',
            color='k', linestyle='dashdot')

    # plt.tight_layout()
    plt.xlabel(r'$E$ (eV)');
    plt.ylabel('Leakage Rate (per source particle per eV)')
    plt.xlim(0, 1.);
    plt.ylim(0, 1.8)
    plt.legend([
        'MCNP (tabular)',
        'MiniMC\n(piecewise constant T)',
        'MiniMC\n(linear T)'], frameon=False)
    plt.grid(True)

    tikzplotlib.clean_figure()
    tikzplotlib.save('spectrum_comparison.tex')

    compare_univariate_pdf(
            f'Beta E={E} eV, T={T} K',
            marginalize(pdf_minimc, axis='beta'),
            marginalize(pdf_mcnp, axis='beta'),
            axis='beta')
    compare_univariate_pdf(
            f'Alpha E={E} eV, T={T} K',
            marginalize(pdf_minimc, axis='alpha'),
            marginalize(pdf_mcnp, axis='alpha'),
            axis='alpha')
