#!/usr/bin/env python
import numpy as np
import pandas as pd
from pyminimc import util
from matplotlib import pyplot as plt
import tikzplotlib

prefix = "/Users/atumulak/Developer/minimc/data/tnsl/endfb8/"
alpha_cutoff = 1636.7475317348378


def quadratic_plot(alpha_df, col_index, p0=None, **kwargs):
    cdfs = alpha_df.index.values
    pdfs = compute_pdf(alpha_df, p0=p0).iloc[:, col_index].values
    alpha_boundaries = alpha_df.iloc[:, col_index].values
    alphas = np.linspace(0, alpha_boundaries[-1] * 0.99, 10000)
    # find index of given alpha value
    a_hi_i = np.searchsorted(alpha_boundaries, alphas, side="right")
    a_hi = alpha_boundaries[a_hi_i]
    a_lo = alpha_boundaries[a_hi_i - 1]
    f_hi = pdfs[a_hi_i]
    f_lo = pdfs[a_hi_i - 1]
    r = (alphas - a_lo) / (a_hi - a_lo)
    F_lo = cdfs[a_hi_i - 1]
    interped = F_lo + (a_hi - a_lo) * (f_lo * r + 0.5 * (f_hi - f_lo) * r * r)
    # print(f"badness: {np.sum(np.square(np.diff(pdfs)))}")
    plt.plot(alphas, interped, label=f"f0 = {pdfs[0]}", **kwargs)
    # plt.scatter(alpha_boundaries, np.concatenate([[0], np.square(np.diff(pdfs))]))
    plt.scatter(alpha_boundaries, cdfs, color="k")


def compute_pdf(alpha, p0=None):
    # compute optimal value of first pdf
    cs = 2 / (alpha.diff().iloc[1:, :].div(np.diff(alpha.index), axis="rows"))
    # number of PDF values that depend on p0 (except p0)
    N = cs.shape[0]

    if p0 is None:
        sum_xi = (
            cs.multiply([(-1) ** i for i in range(N)], axis="rows")
            .cumsum()
            .iloc[:-1, :]
        )
        square_diff_alpha = np.square(alpha.shift(-1) - alpha).iloc[1:-1]
        p0 = (square_diff_alpha * sum_xi).sum() / (square_diff_alpha.sum())

    pdf_df = (
        cs.multiply([(-1) ** (i) for i in range(N)], axis="rows").cumsum() - p0
    ).multiply([(-1) ** i for i in range(N)], axis="rows")
    pdf_df.loc[0.0] = p0
    pdf_df = pdf_df.sort_index()
    return pdf_df


for i in range(1, 4):
    U_df = pd.read_hdf(prefix + f"alpha_endfb8_{i}_CDF_coeffs.hdf5")
    S_df = pd.read_hdf(prefix + f"alpha_endfb8_{i}_S_coeffs.hdf5")
    V_df = pd.read_hdf(prefix + f"alpha_endfb8_{i}_beta_T_coeffs.hdf5")
    alpha_df = util.from_svd_dfs(U_df, S_df, V_df)
    # remove rows which have any negative alphas
    is_all_nonneg_alpha = (alpha_df >= 0).all(axis="columns")
    alpha_df = alpha_df[is_all_nonneg_alpha]
    alpha_df.loc[0.0, :] = 0.0
    alpha_df.loc[1.0, :] = alpha_cutoff
    alpha_df = alpha_df.sort_index()
    pdf_df = compute_pdf(alpha_df)
    while (pdf_df < 0).any().any():
        negative_pdf_count = (pdf_df < 0).any(axis="columns").sum()
        print(
            f"negative PDF count: {negative_pdf_count}, alpha_df.shape: {alpha_df.shape}"
        )
        worst_idx = pdf_df[pdf_df < 0].idxmin().mode()[0]
        if worst_idx == pdf_df.index[0]:
            worst_idx = pdf_df.index[1]
        elif worst_idx == pdf_df.index[-1]:
            worst_idx = pdf_df.index[-2]
        # print(f"dropping {worst_idx}")
        alpha_df = alpha_df.drop(worst_idx)
        pdf_df = compute_pdf(alpha_df)
        # quadratic_plot(alpha_df, col_index=10)
        # quadratic_plot(alpha_df, col_index=20)
        quadratic_plot(alpha_df, col_index=30, color="k")
        quadratic_plot(alpha_df, col_index=30, color="k", p0=0.6)
        plt.xlabel(r"$\alpha$")
        plt.ylabel(r"CDF")
        plt.xlim([0, 25])
        plt.ylim([0, 1.1])
        plt.xlim(left=0)
        # plt.legend()
        tikzplotlib.clean_figure()
        tikzplotlib.save(
            "/Users/atumulak/Developer/phd-dissertation/figures/piecewise-linear.tex"
        )
        exit()
    # while not (pdf_df >= 0).all().all():
    #     # remove first row which has negative PDF
    #     is_all_nonneg_pdf = (pdf_df >= 0).all(axis="columns")
    #     # get first row where a negative PDF is found
    #     first_nonneg_pdf_index = is_all_nonneg_pdf[~is_all_nonneg_pdf].index[0]
    #     alpha_df = alpha_df.drop(first_nonneg_pdf_index)
    #     pdf_df = compute_pdf(alpha_df)
