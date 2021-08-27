#include "ThermalScattering.hpp"
#include "Cell.hpp"
#include "Constants.hpp"
#include "Particle.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <random>
#include <variant>

// ThermalScattering

//// public

ThermalScattering::ThermalScattering(const pugi::xml_node& tsl_node) noexcept
    : beta_cdf{tsl_node.attribute("beta_cdf_file").as_string()},
      alpha_cdf{tsl_node.attribute("alpha_cdf_file").as_string()},
      beta_cutoff{tsl_node.attribute("beta_cutoff").as_double()},
      alpha_cutoff{tsl_node.attribute("alpha_cutoff").as_double()},
      min_temperature{tsl_node.attribute("min_temperature").as_double()},
      max_temperature{tsl_node.attribute("max_temperature").as_double()},
      awr{tsl_node.parent().parent().parent().attribute("awr").as_double()}{}


//// private

Real ThermalScattering::EvaluateBeta(
    const size_t E_index, const size_t cdf_index,
    Temperature T) const noexcept {
  Beta result = 0.;
  for (size_t order = 0; order < n_beta_coefficients; order++) {
    result +=
        beta_cdf.at(E_index, cdf_index, order) * std::pow(T, -(order / 2.));
  }
  return result;
}

Real ThermalScattering::EvaluateAlpha(
    const size_t beta_index, const size_t cdf_index,
    Temperature T) const noexcept {
  Alpha result = 0.;
  for (size_t order = 0; order < n_alpha_coefficients; order++){
    result += alpha_cdf.at(beta_index, cdf_index, order) *
              std::pow(T, -static_cast<double>(order));
  }
  return result;
}

ThermalScattering::Beta ThermalScattering::SampleBeta(
    Particle& p, ContinuousEnergy E, Temperature T) const noexcept {

  // find index of energy value strictly greater than E
  const size_t E_hi_i = std::distance(
      beta_energies.cbegin(),
      std::upper_bound(beta_energies.cbegin(), beta_energies.cend(), E));
  // We require that E is strictly less than cutoff energy
  assert (E_hi_i != beta_energies.size());

  // Determine interpolation factor. If E < E_min, we force the use of the grid
  // at E_min by setting r = 1. This way, index of the sampled energy E_s_i is
  // always valid.
  const Real r =
      E_hi_i != 0
          ? (E - beta_energies.at(E_hi_i - 1)) /
                (beta_energies.at(E_hi_i) - beta_energies.at(E_hi_i - 1))
          : 1;
  const size_t E_s_i =
      r <= std::uniform_real_distribution{}(p.rng) ? E_hi_i - 1 : E_hi_i;
  const ContinuousEnergy E_s = beta_energies.at(E_s_i);

  // sample a CDF value
  const Real F = std::uniform_real_distribution{}(p.rng);
  // find index of CDF value strictly greater than sampled CDF value
  const size_t F_hi_i = std::distance(
      beta_cdfs.cbegin(),
      std::upper_bound(beta_cdfs.cbegin(), beta_cdfs.cend(), F));

  // Evaluate nearest Fs on the sampled E grid.
  const Real F_lo = F_hi_i != 0 ? beta_cdfs.at(F_hi_i - 1) : 0;
  const Real F_hi = F_hi_i != beta_cdfs.size() ? beta_cdfs.at(F_hi_i) : 1;

  const ContinuousEnergy E_s_b_lo =
      F_hi_i != 0 ? EvaluateBeta(E_s_i, F_hi_i - 1, T)
                  : -(E_s * 1e6) / (constants::boltzmann * T);
  const ContinuousEnergy E_s_b_hi =
      F_hi_i != beta_cdfs.size() ? EvaluateBeta(E_s_i, F_hi_i, T) : beta_cutoff;

  // Evalute interpolated value of beta on the sampled E grid (assuming
  // histogram PDF)
  const Beta b_prime =
      E_s_b_lo + (F - F_lo) / (F_hi - F_lo) * (E_s_b_hi - E_s_b_lo);

  // Scale to preserve thresholds at actual incident energy E. The minimum and
  // maximum values of beta are not included in the dataset since their
  // analytical form is known.
  const Beta b_min = -(E * 1e6) / (constants::boltzmann * T);
  const Beta b_max = beta_cutoff;
  const Beta E_s_b_min = -(E_s * 1e6) / (constants::boltzmann * T);
  const Beta E_s_b_max = beta_cutoff;
  return b_min + (b_prime - E_s_b_min) * (b_max - b_min) /
                     (E_s_b_max - E_s_b_min);
}

ThermalScattering::Alpha ThermalScattering::SampleAlpha(
    Particle& p, const Beta& b, ContinuousEnergy E,
    Temperature T) const noexcept {

  // assume S(a,b) = S(a,-b)
  const Beta abs_b = std::abs(b);
  // get sign of beta: https://stackoverflow.com/a/4609795/5101335
  const int_fast8_t sgn_b = (0 < b) - (b < 0);
  // find index of beta value strictly greater than abs_b
  const size_t b_hi_i = std::distance(
      alpha_betas.cbegin(),
      std::upper_bound(alpha_betas.cbegin(), alpha_betas.cend(), abs_b));

  // For any given abs_b, we have two beta grids to choose from:
  // b_lo <= abs_b < b_hi
  // For positive b, if beta_cutoff <= b_hi, we force the use of b_lo
  // For negative b, if -b_hi < - E / (k * T), we force the use of b_lo
  const bool snap_to_lower = (
      (sgn_b == 1 && beta_cutoff <= alpha_betas.at(b_hi_i)) ||
      (sgn_b == -1 &&
       -alpha_betas.at(b_hi_i) < -(E * 1e6) / (constants::boltzmann * T)));
  // Determine interpolation factor and sample a beta grid to use, b_s
  const Real r =
    snap_to_lower ? 0 :
    (abs_b - alpha_betas.at(b_hi_i - 1)
     / (alpha_betas.at(b_hi_i) - alpha_betas.at(b_hi_i - 1)));
  const size_t b_s_i =
      r <= std::uniform_real_distribution{}(p.rng) ? b_hi_i - 1 : b_hi_i;
  Beta b_s = sgn_b * alpha_betas.at(b_s_i);

  // compute minimum and maximum possible value of alpha
  const auto sqrt_E = std::sqrt(E * 1e6);
  const auto b_s_sqrt_E_bkT =
      std::sqrt(E * 1e6 + b_s * constants::boltzmann * T);
  const auto akT = awr * constants::boltzmann * T;
  // minimum and maximum alpha values at sampled beta
  const Alpha b_s_a_min = std::pow(sqrt_E - b_s_sqrt_E_bkT, 2) / akT;
  const Alpha b_s_a_max = std::pow(sqrt_E + b_s_sqrt_E_bkT, 2) / akT;
  assert(b_s_a_max < alpha_cutoff);
  // helper function to find the CDF value that would return `a`. (assuming
  // histogram PDF)
  const auto find_cdf = [this, &b_s_i, &T](const Alpha& a) -> Real {
    // find index of CDF value that evaluates to be strictly greater than `a`
    const size_t F_a_hi_i = std::distance(
        alpha_cdfs.cbegin(),
        std::upper_bound(
            alpha_cdfs.cbegin(), alpha_cdfs.cend(), a,
            [this, &b_s_i, &T](const Alpha& a, const Real& cdf) {
              // https://en.cppreference.com/w/cpp/language/operator_arithmetic
              const size_t cdf_i = &cdf - &alpha_cdfs.front();
              return a < EvaluateAlpha(b_s_i, cdf_i, T);
            }));
    // evaluate nearest Fs
    const Real F_a_lo = F_a_hi_i != 0 ? alpha_cdfs.at(F_a_hi_i - 1) : 0;
    const Real F_a_hi =
        F_a_hi_i != alpha_cdfs.size() ? alpha_cdfs.at(F_a_hi_i) : 1;
    // evaluate nearest alphas
    const Alpha a_lo =
        F_a_hi_i != 0 ? EvaluateAlpha(b_s_i, F_a_hi_i - 1, T) : 0;
    const Alpha a_hi = F_a_hi_i != alpha_cdfs.size()
                           ? EvaluateAlpha(b_s_i, F_a_hi_i, T)
                           : alpha_cutoff;
    // assume histogram PDF
    return F_a_lo + (a - a_lo) * (F_a_hi - F_a_lo) / (a_hi - a_lo);
  };
  // CDFs corresponding to b_s_a_min and b_s_a_max
  const auto F_min = find_cdf(b_s_a_min);
  const auto F_max = find_cdf(b_s_a_max);

  // sample a CDF value which is scaled to return a result in [b_s_a_min,
  // b_s_a_max]
  const Real F =
      F_min + std::uniform_real_distribution{}(p.rng) * (F_max - F_min);
  // find index of CDF value strictly greater than sampled CDF value
  const size_t F_hi_i = std::distance(
      alpha_cdfs.cbegin(),
      std::upper_bound(alpha_cdfs.cbegin(), alpha_cdfs.cend(), F));

  // Evaluate nearest Fs on the sampled beta grid
  const Real F_lo = F_hi_i != 0 ? alpha_cdfs.at(F_hi_i - 1) : 0;
  const Real F_hi = F_hi_i != alpha_cdfs.size() ? alpha_cdfs.at(F_hi_i) : 1;

  // Evaluate nearest alphas on the sampled beta grid.
  const Real b_s_a_lo = F_hi_i != 0 ? EvaluateAlpha(b_s_i, F_hi_i - 1, T) : 0;
  const Real b_s_a_hi = F_hi_i != alpha_cdfs.size()
                            ? EvaluateAlpha(b_s_i, F_hi_i, T)
                            : alpha_cutoff;

  // Evalute interpolated value of alpha on the sampled beta grid (assuming
  // histogram PDF)
  const Beta a_prime =
      b_s_a_lo + (F - F_lo) / (F_hi - F_lo) * (b_s_a_hi - b_s_a_lo);

  // Scale to preserve thresholds at actual beta value b. The minimum and
  // maximum values of alpha have analytical forms which are known.
  const auto b_sqrt_E_bkT = std::sqrt(E * 1e6 + b * constants::boltzmann * T);
  const Alpha b_a_min = std::pow(sqrt_E - b_sqrt_E_bkT, 2) / akT;
  const Alpha b_a_max = std::pow(sqrt_E + b_sqrt_E_bkT, 2) / akT;
  return b_a_min +
         (a_prime - b_s_a_min) * (b_a_max - b_a_min) / (b_s_a_max - b_s_a_min);
}
