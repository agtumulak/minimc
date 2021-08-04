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
from multiprocessing import Pool
from scipy import integrate, interpolate, optimize
from tqdm import tqdm


# Raw S(alpha, beta, T) data
realization = 1
df = pd.read_hdf('/Users/atumulak/Developer/sab-poly-fit/sab_medium_retry.hdf5')
# df = pd.read_hdf('/Users/atumulak/Developer/sab-data-generator/full_order200/df.hdf5')

df = (
        df[df['Realization']==realization]
        .drop('Realization', axis='columns')
        .set_index(['beta', 'alpha', 'T']))
# df = pd.read_hdf('../data/H_in_H2O_sab.hdf5')
df_betas = np.array(sorted(df.index.unique(level='beta')))
df_alphas = np.array(sorted(df.index.unique(level='alpha')))
df_Ts = np.array(sorted(df.index.unique(level='T')))

df_interpolator = interpolate.RegularGridInterpolator(
        (df_betas, df_alphas, df_Ts),
        df.values.reshape(len(df_betas), len(df_alphas), len(df_Ts)),
        bounds_error=False, fill_value=None)

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

# dimensionless energy transfer
b_grid = np.array([
    0, 0.006375, 0.01275, 0.0255, 0.03825, 0.051, 0.06575, 0.0806495, 0.120974,
    0.161299, 0.241949, 0.322598, 0.403248, 0.483897, 0.564547, 0.645197,
    0.725846, 0.806496, 0.887145, 0.967795, 1.04844, 1.12909, 1.20974, 1.29039,
    1.37104, 1.45169, 1.53234, 1.61299, 1.69364, 1.77429, 1.85494, 1.93559,
    2.01624, 2.09689, 2.17754, 2.25819, 2.33884, 2.41949, 2.50014, 2.58079,
    2.6695, 2.76709, 2.87445, 2.9925, 3.12235, 3.2653, 3.42247, 3.59536,
    3.78549, 3.99467, 4.22473, 4.47787, 4.75631, 5.06258, 5.39939, 5.76997,
    6.17766, 6.62607, 7.11924, 7.66181, 8.25862, 8.91511, 9.63722, 10.432,
    11.3051, 12.2668, 13.3243, 14.4867, 15.766, 17.1733, 18.7218, 20.4245,
    22.2976, 24.3572, 26.6234, 29.1165, 31.8586, 34.8759, 38.1936, 41.844,
    45.8583, 50.2749, 55.1331, 60.4771, 66.3554, 72.8215, 79.9338, 90, 100, 110,
    120, 130, 140, 150, 160])

# temperature grid in K
T_grid = np.array([
    300, 310, 320, 330, 340, 350, 360, 370, 380, 390, 400, 410, 420, 430, 440,
    450, 460, 470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590,
    600, 610, 620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740,
    750, 760, 770, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890,
    900, 910, 920, 930, 940, 950, 960, 970, 980, 990, 1000])

# value of CDF
F = np.array([
    0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3,
    0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6,
    0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
    0.925, 0.95, 0.975])


def process_E_T(args):
    E, T = args
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
        S_values = df.xs((beta if beta > 0 else -beta, T), level=('beta', 'T'))['S']
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
        alpha_pdfs[beta] = S_values / alpha_integral
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
    beta, T = args
    S_values = df.xs((beta, T), level=('beta', 'T'))['S']
    E_independent_alpha_cdf = pd.Series(
            np.concatenate((
                [0],
                integrate.cumulative_trapezoid(S_values, S_values.index))),
                index=S_values.index)
    E_independent_alpha_integral = E_independent_alpha_cdf.iloc[-1]
    return (E_independent_alpha_cdf / E_independent_alpha_integral).values


def plot_beta_pdf(E, T):
    beta_pdf_pdos, _, _, _ = process_E_T((E, T))
    beta_alpha_pdf_mcnp = pd.read_hdf('~/Downloads/MCNP6/hockey_stick_E_mu.hdf5', 'hockey_stick_E_mu')
    beta_alpha_pdf_mcnp['absolute error'] = (
            beta_alpha_pdf_mcnp['estimate'] * beta_alpha_pdf_mcnp['relative error'])
    beta_pdf_mcnp = beta_alpha_pdf_mcnp.groupby(level='energy').aggregate({
        'estimate': sum,
        'absolute error': lambda s: sum(s ** 2) ** 0.5 }) # assume uncorrelated
    beta_pdf_mcnp.index = beta_pdf_mcnp.index * 1e6 # MeV to eV
    # convert mcnp estimates to probability density in beta.
    beta_pdf_mcnp['bin width'] = (
            beta_pdf_mcnp.index - np.concatenate(([0], beta_pdf_mcnp.index[:-1])))
    beta_pdf_mcnp['density'] = beta_pdf_mcnp['estimate'] / beta_pdf_mcnp['bin width'] * (k * T)
    beta_pdf_mcnp['density error'] = beta_pdf_mcnp['absolute error'] / beta_pdf_mcnp['bin width'] * (k * T)
    beta_pdf_mcnp = pd.concat((
        pd.DataFrame(0, index=[0.], columns=beta_pdf_mcnp.columns),
        beta_pdf_mcnp))
    # Assumes beta boundaries correspond to energy bins in mcnp.
    df = pd.DataFrame({
        'mcnp density': beta_pdf_mcnp['density'].values,
        'mcnp density lower': beta_pdf_mcnp['density'].values - beta_pdf_mcnp['density error'].values,
        'mcnp density upper': beta_pdf_mcnp['density'].values + beta_pdf_mcnp['density error'].values,
        'pdos density': beta_pdf_pdos.values},
        index=beta_pdf_pdos.index)
    plt.plot(df.index, df['mcnp density'], label='mcnp')
    plt.fill_between(df.index, df['mcnp density lower'], df['mcnp density upper'], alpha=0.25, color='C0')
    plt.plot(df.index, df['pdos density'], label='pdos')
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$p_{\beta}(\beta)$')
    plt.xlim(df.index.min(), df.index.max())
    plt.ylim(0, df[['mcnp density', 'pdos density']].max().max())
    plt.legend()
    plt.show()


def get_alpha_pdf(df, b_lower=None, b_upper=None):
    '''
    Columns are beta values. Rows are alpha values. All entries are valid (no
    nans).
    '''
    b_lower = b_lower if b_lower else df.columns.min()
    b_upper = b_upper if b_upper else df.columns.max()
    # add columns for lower and upper beta
    if b_lower not in df.columns:
        df[b_lower] = np.nan
    if b_upper not in df.columns:
        df[b_upper] = np.nan
    alpha_pdf = (
            df[sorted(df.columns)]
            .interpolate(method='index', axis='columns')
            .loc[:, b_lower:b_upper]
            .apply(
                lambda s:
                integrate.trapezoid(s, s.index),
                axis='columns'))
    alpha_integral = integrate.trapezoid(alpha_pdf, alpha_pdf.index)
    alpha_pdf /= alpha_integral
    return alpha_pdf


def get_alpha_pdf_mcnp(E, T, b_lower=None, b_upper=None):
    df = pd.read_hdf('~/Downloads/MCNP6/hockey_stick_E_mu.hdf5', 'hockey_stick_E_mu')
    df['absolute error'] = df['estimate'] * df['relative error']
    df.index = df.index.set_levels(
            df.index.unique(level='energy') * 1e6, level='energy')
    # compute bin widths
    cosines = np.array(sorted(df.index.unique(level='cosine')))
    cosine_widths = cosines - np.concatenate(([-1.0], cosines[:-1]))
    Es = np.array(sorted(df.index.unique(level='energy')))
    E_widths = Es - np.concatenate(([0.0], Es[:-1]))
    # compute probability densities
    df['cosine'] = df.index.unique(level='cosine')
    df['energy'] = df.index.unique(level='energy')
    df['cosine widths'] = df['cosine']
    df['energy widths'] = df['energy']
    df = df.replace(to_replace={
        'cosine widths': {k: v for k, v in zip(cosines, cosine_widths)},
        'energy widths': {k: v for k, v in zip(Es, E_widths)}})
    jacobian = A * (k * T) ** 2 / (2 * np.sqrt(E * df['energy']))
    df['density'] = df['estimate'] / (df['cosine widths'] * df['energy widths']) * jacobian
    df['density error'] = df['absolute error'] / (df['cosine widths'] * df['energy widths']) * jacobian
    # compute corresponding beta values
    df['beta'] = (df['energy'] - E) / (k * T)
    # compute corresponding alpha values. use lower cosine boundaries
    df['lower cosine'] = df['cosine']
    df = df.replace(to_replace={
        'lower cosine': {k: v for k, v in zip(cosines, np.concatenate(([-1.], cosines[:-1])))}})
    df['alpha'] = (E + df['energy'] - 2. * df['cosine'] * (df['energy'] * E) ** 0.5) / (A * k * T)
    b_lower = b_lower if b_lower is not None else df['beta'].min()
    b_upper = b_upper if b_upper is not None else df['beta'].max()
    density_df = (
            df
            .pivot(index='alpha', columns='beta', values='density')
            .interpolate(method='index', axis='index', limit_area='inside')
            .fillna(0))
    alpha_pdf = get_alpha_pdf(density_df, b_lower=b_lower, b_upper=b_upper)
    return alpha_pdf
    # alpha_pdf_error = (
    #         df
    #         .pivot(index='alpha', columns='beta', values='density error')
    #         .loc[:, b_lower:b_upper]
    #         .apply(
    #             lambda row: ((row.dropna() ** 2).sum()) ** 0.5, axis='columns'))
    # plt.plot(alpha_pdf.index, alpha_pdf)
    # plt.fill_between(
    #         alpha_pdf.index,
    #         alpha_pdf - alpha_pdf_error,
    #         alpha_pdf + alpha_pdf_error,
    #         alpha=0.25, color='C0')
    # plt.show()


def get_alpha_pdf_pdos(beta_pdf_pdos, alpha_pdfs_pdos, b_lower=None, b_upper=None):
    beta_alpha_pdf_pdos = (
            pd.DataFrame(alpha_pdfs_pdos, )
            .mul(beta_pdf_pdos.iloc[1:-1])
            .interpolate(method='index', axis='index', limit_area='inside')
            .fillna(0))
    beta_alpha_pdf_pdos.index.rename('alpha')
    return get_alpha_pdf(beta_alpha_pdf_pdos, b_lower=b_lower, b_upper=b_upper)


def compare_alpha_pdf(E, T, b_lower=None, b_upper=None):
    alpha_pdf_pdos = get_alpha_pdf_pdos(
            *process_E_T((E, T))[:2], b_lower=b_lower, b_upper=b_upper)
    alpha_pdf_mcnp = get_alpha_pdf_mcnp(
            E, T, b_lower=b_lower, b_upper=b_upper)
    plt.plot(alpha_pdf_pdos.index, alpha_pdf_pdos, label='pdos')
    plt.plot(alpha_pdf_mcnp.index, alpha_pdf_mcnp, label='mcnp')
    plt.xlabel(r'$\alpha$')
    plt.ylabel(r'$p(\alpha \mid \beta_{\min} < \beta < \beta_{\max})$')
    b_lower = b_lower if b_lower is not None else - E / (k * T)
    b_upper = b_upper if b_upper is not None else 20
    plt.title(r'$\beta_{\min} = ' + str(b_lower) + r', \beta_{\max} = ' + str(b_upper) + r'$')
    plt.legend()
    plt.show()


def vary_T():
    max_cdf = 0
    # for T in df_Ts:
    for T in [273.6]:
        beta, beta_cdf, _, _ = process_E_T((1.0, T))
        max_cdf = max(max_cdf, beta_cdf.max())
        plt.plot(beta, beta_cdf, label=f'{T} K')
    plt.ylim(0, max_cdf)
    plt.xlabel(r'$\beta$')
    plt.ylabel(r'$P_{\beta}(\beta)$')
    plt.legend()
    plt.show()


def beta_functional_expansion(x):
    beta_df = pd.DataFrame(
            0,
            index=pd.Index(F, name='CDF'),
            columns=pd.MultiIndex.from_product((E_grid, df_Ts), names=('E', 'T')))
    for E, x_E in zip(E_grid, x):
        for T, beta_cdf in zip(df_Ts, x_E):
            beta_df.loc[:, (E, T)] = np.interp(F, beta_cdf, beta_cdf.index)
    def fitting_function(T, c0, c1, c2):
        # clamp T to be within [min_val, max_val]
        min_val = 0.5
        max_val = 1.0
        x = (T - T.min()) / (T.max() - T.min()) * (max_val - min_val) + min_val
        return c0 + c1 / np.sqrt(x) + c2 / x
    beta_df_fit = beta_df.unstack().groupby(level=['E', 'CDF']).apply(
            lambda s: pd.Series(
                coeffs := optimize.curve_fit(fitting_function, s.index.unique(level='T'), s)[0],
                index=pd.Index(range(len(coeffs)), name='coefficient')))
    # check that CDFS are monotonic for certain T values
    test_T = np.linspace(df_Ts.min(), df_Ts.max(), 100)
    beta_df_reconstructed = beta_df_fit.groupby(level=['E', 'CDF']).apply(
            lambda s: pd.Series(
                fitting_function(test_T, *s), index=pd.Index(test_T, name='T')))
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


def alpha_functional_expansion(x):
    max_beta = max(20, E_grid[-1] / (k * df_Ts[0]))
    betas = df_betas[:np.searchsorted(df_betas, max_beta)] # will never include max_beta
    alpha_df = pd.DataFrame(
            0,
            index=pd.Index(F, name='CDF'),
            columns=pd.MultiIndex.from_product((betas, df_Ts), names=('beta', 'T')))
    for beta, x_beta in zip(betas, x):
        for T, alpha_cdf in zip(df_Ts, x_beta):
            alpha_df.loc[:, (beta, T)] = np.interp(F, alpha_cdf, df_alphas)
    def fitting_function(T, c0, c1, c2):
        # clamp T to be within [min_val, max_val]
        min_val = 0.5
        max_val = 1.0
        x = (T - T.min()) / (T.max() - T.min()) * (max_val - min_val) + min_val
        return c0 + c1 / x + c2 / x ** 2
    alpha_df_fit = alpha_df.unstack().groupby(level=['beta', 'CDF']).apply(
            lambda s: pd.Series(
                coeffs := optimize.curve_fit(fitting_function, s.index.unique(level='T'), s)[0],
                index=pd.Index(range(len(coeffs)), name='coefficient')))
    # check that CDFS are monotonic for certain T values
    test_T = np.linspace(df_Ts.min(), df_Ts.max(), 100)
    alpha_df_reconstructed = alpha_df_fit.groupby(level=['beta', 'CDF']).apply(
            lambda s: pd.Series(
                fitting_function(test_T, *s), index=pd.Index(test_T, name='T')))
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


def process_all_E_T():
    with Pool(processes=10) as pool:
        # incident energy, temperature pairs
        E_T_values = np.array([(E, T) for E in E_grid for T in df_Ts])
        return (
                np.array(
                    [x[2] for x in tqdm(
                        pool.imap(func=process_E_T, iterable=E_T_values),
                        total=len(E_T_values))],
                    dtype=object).reshape(len(E_grid), len(df_Ts)))


def process_all_b_T():
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
