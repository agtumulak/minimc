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
from multiprocessing import Pool
from scipy import integrate, interpolate, optimize
from tqdm import tqdm


# target nuclide mass ratio
A = 0.999167339

# Boltzmann constant in eV / K
k = 8.617333262145E-5

# incident neutron energy grid in eV (for H in H2O)
E_grid = np.array([
    1.00E-5, 1.78E-5, 2.50E-5, 3.50E-5, 5.0E-5, 7.00E-5, 1.00E-4, 1.26E-4,
    2.00E-4, 2.53E-4, 2.97E-4, 3.50E-4, 4.20E-4, 5.06E-4, 6.15E-5, 7.50E-4,
    8.70E-4, 1.01E-3, 1.23E-3, 1.50E-3, 1.80E-3, 2.03E-3, 2.28E-3, 2.60E-3,
    3.00E-3, 3.50E-3, 4.05E-3, 4.50E-3, 5.00E-3, 5.60E-3, 6.33E-3, 7.20E-3,
    8.10E-3, 9.11E-3, 1.00E-2, 1.06E-2, 1.15E-2, 1.24E-2, 1.33E-2, 1.42E-2,
    1.50E-2, 1.62E-2, 1.82E-2, 1.99E-2, 2.05E-2, 2.15E-2, 2.28E-2, 2.53E-2,
    2.80E-2, 3.06E-2, 3.38E-2, 3.65E-2, 3.95E-2, 4.28E-2, 4.65E-2, 5.00E-2,
    5.69E-2, 6.25E-2, 6.90E-2, 7.50E-2, 8.20E-2, 9.00E-2, 9.60E-2, 1.04E-1,
    1.12E-1, 1.20E-1, 1.28E-1, 1.36E-1, 1.46E-1, 1.60E-1, 1.72E-1, 1.84E-1,
    2.00E-1, 2.28E-1, 2.51E-1, 2.71E-1, 2.91E-1, 3.01E-1, 3.21E-1, 3.58E-1,
    3.90E-1, 4.17E-1, 4.50E-1, 5.03E-1, 5.60E-1, 6.25E-1, 7.00E-1, 7.80E-1,
    8.60E-1, 9.50E-1, 1.05, 1.16, 1.28, 1.42, 1.55])

# value of CDF
F = np.array([
    0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3,
    0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6,
    0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
    0.925, 0.95, 0.975])


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
            if sab_file.readline().endswith('1 7  4    5\n'):
                N_beta = int(sab_file.readline().split()[0])
                break
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
        return (
                pd.DataFrame.from_dict(
                    {'alpha': alphas, 'beta': betas, 'T': Ts, 'S': Ss})
                .set_index(['beta', 'alpha', 'T']))


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
    df_betas = np.array(sorted(sab_df.index.unique('beta')))
    # valid beta values
    min_beta = - E / (k * T)
    max_beta = 20
    min_beta_index = np.searchsorted(df_betas, -min_beta) # as range end, will never include min_beta
    max_beta_index = np.searchsorted(df_betas, max_beta) # as range end, will never include max_beta
    betas = np.concatenate((
        [min_beta],
        -np.flip(df_betas[1:min_beta_index]),
        df_betas[:max_beta_index],
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
        alpha_cdf = pd.Series(
                np.concatenate((
                    [0],
                    integrate.cumulative_trapezoid(S_values, S_values.index))),
                index=S_values.index)
        alpha_integral = alpha_cdf.iloc[-1]
        alpha_cdfs[beta] = alpha_cdf / alpha_integral
        alpha_pdfs[beta] = (S_values / alpha_integral).rename(beta)
        beta_pdf.loc[beta] = np.exp(-beta / 2.0) * alpha_integral
    # convert pdf to cdf
    beta_cdf = pd.Series(
            np.concatenate((
                [0],
                integrate.cumulative_trapezoid(beta_pdf, beta_pdf.index))),
            index=beta_pdf.index)
    beta_integral = beta_cdf.iloc[-1]
    beta_cdf /= beta_integral
    beta_pdf /= beta_integral
    return beta_pdf, alpha_pdfs, beta_cdf, alpha_cdfs


def process_b_T(args):
    sab_df, beta, T = args
    S_values = sab_df.xs((beta, T), level=('beta', 'T'))['S']
    E_independent_alpha_cdf = pd.Series(
            np.concatenate((
                [0],
                integrate.cumulative_trapezoid(S_values, S_values.index))),
                index=S_values.index)
    E_independent_alpha_integral = E_independent_alpha_cdf.iloc[-1]
    return (E_independent_alpha_cdf / E_independent_alpha_integral).values


def beta_fitting_function(T, c0, c1, c2, c3, c4, c5):
    return c0 + c1 / T ** (1/2) + c2 / T + c3 / T ** (3/2) + c4 / T ** 2 + c5 / T ** (5/2)


def alpha_fitting_function(T, c0, c1, c2):
    return c0 + c1 / T + c2 / T ** 2


def get_pdf_pdos(sab_df, E, T):
    """
    Creates bivariate PDF in alpha and beta from given DataFrame

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
    beta_pdf, alpha_pdfs, _, _ = process_E_T((sab_df, E, T))
    for beta, p_beta in beta_pdf.iloc[1:-1].iteritems():
        alpha_pdfs[beta] *= p_beta
    return (
            pd
            .concat(alpha_pdfs.values(), names=alpha_pdfs.keys(), axis='columns')
            .interpolate(method='index', axis='index', limit_area='inside')
            .fillna(0)
            .stack()
            .rename_axis(['alpha', 'beta'])
            .rename('pdos'))


def get_pdf_pdos_reconstructed(beta_hdf_path, alpha_hdf_path, E, T):
    """
    Reconstructs bivariate PDF in alpha and beta from expansion
    coefficients in temperature for alpha and beta

    Parameters
    ---------
    beta_hdf_path : string
        Path to CDF file for beta
    alpha_cdf_path : string
        Path to conditional CDF file for alpha
    E : float
        Incident energy in eV
    T : float
        Temperature in K
    """
    # evaluate beta at E and T
    beta_cdf = (
            pd.read_hdf(beta_hdf_path)
            .loc[E * 1e-6]
            .groupby('CDF')
            .apply(lambda s: beta_fitting_function(T, *s.values)[0]))
    # add endpoints at F = 0 and F = 1
    beta_cdf = pd.concat((
        pd.Series({0: - E / (k * T)}),
        beta_cdf,
        pd.Series({1: 20})))
    # compute PDF
    beta_pdf = pd.Series(
            (beta_cdf.index[1:] - beta_cdf.index[:-1])
          / (beta_cdf.values[1:] - beta_cdf.values[:-1]),
            index=beta_cdf.values[:-1])
    # add endpoints
    beta_pdf = pd.concat((
        pd.Series({- E / (k * T): 0}),
        beta_pdf,
        pd.Series({20: 0})))
    # evaluate alpha at E and k
    alpha_cdf = (
            pd
            .read_hdf(alpha_hdf_path)
            .groupby(['beta', 'CDF'])
            .apply(lambda s: alpha_fitting_function(T, *s.values)[0]))
    # append alpha CDFs for negative beta values
    alpha_betas = alpha_cdf.index.unique('beta')
    # find largest beta in alpha_betas which is strictly less than E / (k * T)
    # we assume beta = 0 exists so result of searchsorted is >= 1
    min_beta = alpha_betas[np.searchsorted(alpha_betas, -beta_cdf.iloc[0])] - 1
    neg_b_alpha_cdf = (
            alpha_cdf
            .loc[alpha_betas[1]:min_beta] # don't include beta = 0
            .rename(index=lambda x: -x, level='beta') # make beta labels negative
            .sort_index(level='beta'))
    # find largest beta in alpha_betas which is strictly less than 20
    max_beta = alpha_betas[np.searchsorted(alpha_betas, 20)] - 1
    alpha_cdf = pd.concat((neg_b_alpha_cdf, alpha_cdf.loc[:max_beta]))
    def get_joint_probability(s):
        """
        Multiplies conditional probability in alpha with probability in beta
        """
        # choose correct subset of CDF values within min_alpha and max_alpha
        nonlocal E
        beta = s.name
        s = pd.Series(s.index.get_level_values('CDF'), index=pd.Index(s, name='alpha'))
        # insert values for min_alpha and max_apha
        min_alpha = np.square(np.sqrt(E) - np.sqrt(E + beta * k * T)) / (A * k * T)
        max_alpha = np.square(np.sqrt(E) + np.sqrt(E + beta * k * T)) / (A * k * T)
        s.loc[min_alpha] = np.nan # interpolate value later
        s.loc[max_alpha] = np.nan # interpolate value later
        s.loc[0] = 0
        s.loc[632.9] = 1
        s = s.sort_index().interpolate(method='index').loc[min_alpha:max_alpha]
        # rescale CDF to be 0 at min_alpha and 1 at max_alpha
        s = (s - s.min()) / (s.max() - s.min())
        alpha_pdf = pd.Series(
                (s.values[1:] - s.values[:-1]) / (s.index[1:] - s.index[:-1]),
                index=s.index[:-1])
        alpha_pdf[min_alpha] = 0
        alpha_pdf[max_alpha] = 0
        alpha_pdf = alpha_pdf.sort_index()
        # interpolate value of beta pdf
        nonlocal beta_pdf
        beta_pdf_value = np.interp(beta, beta_pdf.index, beta_pdf, left=0, right=0)
        return alpha_pdf * beta_pdf_value
    return (
            alpha_cdf
            .groupby('beta')
            .apply(get_joint_probability)
            .unstack('beta')
            .interpolate(method='index', axis='index', limit_area='inside')
            .fillna(0)
            .stack()
            .rename_axis(['alpha', 'beta'])
            .rename('pdos reconstructed'))


def get_pdf_minimc(counts_path, *bounds_paths):
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


def get_pdf_mcnp(mctal_path, E, T):
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
    # convert to alpha, beta space
    def myfoo(s):
        beta = s.name
        out_E = E + beta * k * T
        mus = s.index.get_level_values('mu')
        alphas = (E + out_E - 2 * mus * np.sqrt(E * out_E)) / (A * k * T)
        return pd.Series(s.values, index=alphas).rename_axis(index={'mu': 'alpha'})
    s = pd.Series(
            density,
            index=pd.MultiIndex.from_product(
                [cosine_bounds[:-1], energy_bounds[:-1]], names=['mu', 'E']),
            name='mcnp')
    # multiply by jacobian
    s = s * 0.5 * A * (k * T) ** 2 / (np.sqrt(s.index.get_level_values('E') * E))
    s = (
            s
            .rename(index=lambda out_E: (out_E - E) / (k * T), level='E')
            .rename_axis(index={'E': 'beta'})
            .groupby('beta')
            .apply(myfoo))
    return (
            s[~s.index.duplicated()]
            .unstack('beta')
            .interpolate(method='index', axis='index', limit_area='inside')
            .fillna(0)
            .stack()
            .rename_axis(['alpha', 'beta'])
            .rename('mcnp'))


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
            if (j == 1):
                ax.set_xlabel(r'$\beta$')
            if (i == 0):
                ax.set_ylabel(r'$\alpha$')
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


def parallel_apply(df_grouped, func):
    with Pool(10) as p:
        results = p.map(func, [(name, group) for name, group in df_grouped], chunksize=100)
    return pd.concat(results)


def fit_points(args):
    group_name, s = args
    s = pd.Series(
            coeffs := optimize.curve_fit(
                beta_fitting_function, s.index.unique(level='T'), s)[0],
            index=pd.Index(range(len(coeffs)), name='coefficient'))
    return pd.concat([s], keys=[group_name], names=['E', 'CDF'])


def beta_functional_expansion(sab_df, x):
    df_Ts = np.array(sorted(sab_df.index.unique('T')))
    beta_df = pd.DataFrame(
            np.nan,
            index=pd.Index(F, name='CDF'),
            columns=pd.MultiIndex.from_product((E_grid, df_Ts), names=('E', 'T')))
    for E, x_E in zip(E_grid, x):
        for T, beta_cdf in zip(df_Ts, x_E):
            beta_df.loc[:, (E, T)] = np.interp(F, beta_cdf, beta_cdf.index)
    beta_df_fit = parallel_apply(
            beta_df.unstack().groupby(level=['E', 'CDF']), fit_points)
    # check that CDFS are monotonic for certain T values
    test_T = np.linspace(df_Ts.min(), df_Ts.max(), 100)
    beta_df_reconstructed = beta_df_fit.groupby(level=['E', 'CDF']).apply(
            lambda s: pd.Series(
                beta_fitting_function(test_T, *s), index=pd.Index(test_T, name='T')))
    print(f"RMSE: {np.sqrt(((beta_df - beta_df_reconstructed.unstack(level='E').unstack(level='T')) ** 2).mean().mean())}")
    is_monotonic = (
            beta_df_reconstructed
            .groupby(level=['E', 'T'])
            .apply(lambda s: s.is_monotonic))
    if not is_monotonic.all():
        print("The following CDFs are not monotonic:")
        print((
            beta_df_reconstructed
            .reorder_levels(['E', 'T', 'CDF'])[~is_monotonic]
            .unstack('E').unstack('T')))
    # set energy units to MeV
    beta_df_fit.index = (
            beta_df_fit.index
            .set_levels(beta_df_fit.index.unique(level='E') * 1e-6, level='E'))
    # return minimc-style array
    return beta_df_fit.to_frame('coefficients').sort_index(level=['E', 'CDF', 'coefficient'])


def alpha_functional_expansion(sab_df, x):
    df_betas = np.array(sorted(sab_df.index.unique('beta')))
    df_alphas = np.array(sorted(sab_df.index.unique('alpha')))
    df_Ts = np.array(sorted(sab_df.index.unique('T')))
    max_beta = max(20, E_grid[-1] / (k * df_Ts[0]))
    betas = df_betas[:np.searchsorted(df_betas, max_beta)] # will never include max_beta
    alpha_df = pd.DataFrame(
            0,
            index=pd.Index(F, name='CDF'),
            columns=pd.MultiIndex.from_product((betas, df_Ts), names=('beta', 'T')))
    for beta, x_beta in zip(betas, x):
        for T, alpha_cdf in zip(df_Ts, x_beta):
            alpha_df.loc[:, (beta, T)] = np.interp(F, alpha_cdf, df_alphas)
    alpha_df_fit = alpha_df.unstack().groupby(level=['beta', 'CDF']).apply(
            lambda s: pd.Series(
                coeffs := optimize.curve_fit(alpha_fitting_function, s.index.unique(level='T'), s)[0],
                index=pd.Index(range(len(coeffs)), name='coefficient')))
    # check that CDFS are monotonic for certain T values
    test_T = np.linspace(df_Ts.min(), df_Ts.max(), 100)
    alpha_df_reconstructed = alpha_df_fit.groupby(level=['beta', 'CDF']).apply(
            lambda s: pd.Series(
                alpha_fitting_function(test_T, *s), index=pd.Index(test_T, name='T')))
    print(f"RMSE: {np.sqrt(((alpha_df - alpha_df_reconstructed.unstack(level='beta').unstack(level='T')) ** 2).mean().mean())}")
    is_monotonic = (
            alpha_df_reconstructed
            .groupby(level=['beta', 'T'])
            .apply(lambda s: s.is_monotonic)
            .all())
    if not is_monotonic:
        print("The following CDFs are not monotonic:")
        print((
            alpha_df_reconstructed
            .reorder_levels(['beta', 'T', 'CDF'])[~is_monotonic]
            .unstack('beta').unstack('T')))
    return alpha_df_fit.to_frame('coefficients').sort_index(level=['beta', 'CDF', 'coefficient'])


def process_all_E_T(sab_df):
    df_Ts = np.array(sorted(sab_df.index.unique('T')))
    with Pool(processes=10) as pool:
        # incident energy, temperature pairs
        E_T_values = np.array([(sab_df, E, T) for E in E_grid for T in df_Ts])
        return (
                np.array(
                    [x[2] for x in tqdm(
                        pool.imap(func=process_E_T, iterable=E_T_values),
                        total=len(E_T_values))],
                    dtype=object).reshape(len(E_grid), len(df_Ts)))


def process_all_b_T(sab_df):
    df_betas = np.array(sorted(sab_df.index.unique('beta')))
    with Pool(processes=10) as pool:
        # largest possible beta value
        max_beta = max(20, E_grid[-1] / (k * df_Ts[0]))
        betas = df_betas[:np.searchsorted(df_betas, max_beta)] # will never include max_beta
        b_T_values = np.array([(b, T) for b in betas for T in df_Ts])
        return (
                np.vstack([
                    x for x in tqdm(
                        pool.imap(func=process_b_T, iterable=b_T_values),
                        total=len(b_T_values))])
                    .reshape(len(betas), len(df_Ts), -1))
