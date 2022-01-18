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
from functools import partial
from inspect import signature
from multiprocessing import Pool
from scipy import optimize
from tqdm import tqdm


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


def process_b_T(args):
    """
    Generates conditional CDF in alpha given beta at a particular value of
    beta and temperature. Returns None if there is zero probability of the
    given beta being sampled.

    Parameters
    ----------
    args: tuple of (group key, group, double)
        The first element is the group key (a (beta, temperature) pair).
        The second element is a pd.Series containing a MultiIndex. The
        MultiIndex levels are temperature `T`, `beta`, and `alpha`. There is
        only a single value of `T` and `beta` while multiple `alpha` values
        must be present. The values are corresponding value of S(a,b,T).
        The third element is the maximum alpha in the entire data set.
    """
    (beta, T), sab_s, max_alpha = args
    # set endpoints to zero
    if not sab_s.index.isin([(T, beta, 0)]).any():
        sab_s[T, beta, 0] = 0
    if not sab_s.index.isin([(T, beta, max_alpha)]).any():
        sab_s[T, beta, max_alpha] = 0
    sab_s = sab_s.sort_index()
    E_independent_alpha_cdf = lin_log_cum_trapz(sab_s[T, beta])
    E_independent_alpha_integral = E_independent_alpha_cdf.iloc[-1]
    # if integral is zero, return nothing
    if E_independent_alpha_integral == 0:
        E_independent_alpha_cdf = None
    else:
        E_independent_alpha_cdf /= E_independent_alpha_integral
        E_independent_alpha_cdf.index = pd.MultiIndex.from_product(
                [[beta], [T], E_independent_alpha_cdf.index],
                names=['beta', 'T', 'alpha'])
    return E_independent_alpha_cdf


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
    # compute PDF
    beta_pdf = pd.Series(
            (beta_cdf.index[1:] - beta_cdf.index[:-1])
          / (beta_cdf.values[1:] - beta_cdf.values[:-1]),
            index=beta_cdf.values[:-1])
    # evaluate alpha at E and T
    alpha_cdf = (
            pd
            .read_hdf(alpha_hdf_path)
            .groupby(['beta', 'CDF'])
            .apply(lambda s: alpha_fitting_function(T, *s.values)[0]))
    # append alpha CDFs for negative beta values
    alpha_betas = alpha_cdf.index.unique('beta')
    # find largest beta in alpha_betas which is strictly less than E / (k * T)
    # we assume beta = 0 exists so result of searchsorted is >= 1
    min_beta = alpha_betas[np.searchsorted(alpha_betas, -beta_cdf.iloc[0]) - 1]
    neg_b_alpha_cdf = (
            alpha_cdf
            .loc[alpha_betas[1]:min_beta] # don't include beta = 0
            .rename(index=lambda x: -x, level='beta') # make beta labels negative
            .sort_index(level='beta'))
    # find largest beta in alpha_betas which is strictly less than 20
    max_beta = alpha_betas[np.searchsorted(alpha_betas, 20) - 1]
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
        # TODO: Get maximum rescaled alpha value programatically
        s.loc[632.9 * 293.6 / 273.6] = 1 # 273.6 K is smallest T in dataset
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


def beta_functional_expansion(sab_df, E_min=1e-5, E_max=4.0, n_Es=100,
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
        Expansion order for proper orthogonal decomposition

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
    # take the union of all CDF values that appear across all incident energies
    all_cdfs = sorted(set(np.concatenate(beta_cdfs.reshape(-1))))
    # choose approximate number of CDF points we want to use
    F = all_cdfs[::len(all_cdfs) // n_cdfs]
    if F[-1] != 1:
        F.append(1)
    print(f"using {len(F)} CDF values")
    beta_df = pd.DataFrame(
            np.nan,
            index=pd.Index(F, name='CDF'),
            columns=pd.MultiIndex.from_product((Es, df_Ts), names=('E', 'T')))
    for E, x_E in zip(Es, beta_cdfs):
        for T, beta_cdf in zip(df_Ts, x_E):
            beta_df.loc[:, (E, T)] = np.interp(F, beta_cdf, beta_cdf.index)
    # perform proper orthogonal decomposition
    beta_df_pod_form = beta_df.stack('T').unstack('CDF')
    U, S, Vt = np.linalg.svd(beta_df_pod_form, full_matrices=False)
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


def alpha_functional_expansion(sab_df, n_betas=100, n_cdfs=1000, order=3):
    """
    Computes the energy-independent conditional CDF in alpha given beta at
    various beta values and temperatures, then performs a functional expansion
    in temperature at various beta values and CDF values.

    Parameters
    ----------
    sab_df : pd.DataFrame
        S(a,b,T) DataFrame
    n_betas : int, optional
        Approximate number of betas to use
    n_cdfs : int, optional
        Approximate number of CDF values to use
    order : int, optional
        Expansion order for proper orthogonal decomposition
    """
    df_Ts = np.array(sorted(sab_df.index.unique('T')))

    # take the union of all beta values that appear across all temperatures
    all_betas = sab_df.index.unique('beta')
    # choose number of beta points we want to use
    selected_betas = all_betas[::len(all_betas) // n_betas]
    # add largest beta from each temperature
    selected_betas = (
            pd.Index(sab_df.groupby('T').apply(
                lambda s: s.index.unique('beta').max()),
                name='beta')
            .union(selected_betas))
    print(f"using {len(selected_betas)} beta values...")
    def common_beta_grid(group):
        """
        Modifies beta values at this temperature to conform with common beta
        grid.
        """
        nonlocal sab_df
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
    alpha_cdfs = (
            parallel_apply(
                common_beta_sab_df.groupby(['beta', 'T']),
                process_b_T, largest_alpha)
            .sort_index())
    # take the union of all CDF values that appear across all incident energies
    all_cdfs = sorted(set(alpha_cdfs))
    # choose number of CDF points we want to use
    F = all_cdfs[::len(all_cdfs) // n_cdfs]
    # don't include 0
    F = F if F[0] != 0 else F[1:]
    # last CDF must always be 1.
    print(f"using {len(F)} CDF values...")
    # interpolate alpha values at selected CDF values
    alpha_df = alpha_cdfs.groupby(['beta', 'T']).apply(
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
        A = U[:, :order] @ np.diag(S[:order]) @ Vt[:order,:]
        alpha_df_pod_form = alpha_df_pod_form.where(~nan_entries, A)
        trace_norm = S.sum()
        abs_rel_diff = np.abs((trace_norm - prev_trace_norm) / prev_trace_norm)
        print (f"trace norm abs rel diff: {abs_rel_diff}", end=" " *  10 + '\r')
        if abs_rel_diff <  1e-7:
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
            columns=alpha_df_pod_form.columns)[~nan_entries]
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


if __name__ == '__main__':
    E = 0.56 # energy in eV
    T = 293.6 # temperature in K
    pdf_pdos = get_pdf_pdos(
            parse_file7('/Users/atumulak/Developer/sab-data-generator/chapman_data.mf7'), E, T)
    pdf_pdos_reconstructed = get_pdf_pdos_reconstructed(
            './beta_cdf_chapman_data.hdf5',
            './alpha_cdf_chapman_data.hdf5',
            E, T)
    pdf_minimc = get_pdf_minimc(
            '../chapman_0.56e-6eV_293.6K.out',
            '../alpha_boundaries.txt',
            '../beta_boundaries.txt')
    pdf_mcnp = get_pdf_mcnp('/Users/atumulak/Downloads/MCNP6/endf80-0.56e-6eV-293.6K.mctal', E, T, label='mcnp (endf80)')
    compare_univariate_pdf(
            f'E={E} eV, T={T} K',
            marginalize(pdf_pdos, axis='beta'),
            marginalize(pdf_pdos_reconstructed, axis='beta'),
            marginalize(pdf_minimc, axis='beta'),
            marginalize(pdf_mcnp, axis='beta'),
            axis='beta')
    compare_univariate_pdf(
            f'E={E} eV, T={T} K',
            marginalize(pdf_pdos, axis='alpha'),
            marginalize(pdf_pdos_reconstructed, axis='alpha'),
            marginalize(pdf_minimc, axis='alpha'),
            marginalize(pdf_mcnp, axis='alpha'),
            axis='alpha')
    compare_bivariate_pdf(
            f'E={E} eV, T={T} K',
            pdf_pdos,
            pdf_pdos_reconstructed,
            pdf_minimc,
            pdf_mcnp)
