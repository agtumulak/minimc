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
  for (size_t order = 0; order <= n_alpha_coefficients; order++){
    result += alpha_cdf.at(beta_index, cdf_index, order) * std::pow(T, -order);
  }
  return result;
}

ThermalScattering::Beta
ThermalScattering::SampleBeta(Particle& p) const noexcept {

  // incident energy of neutron
  const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
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

  // Evaluate nearest betas on the sampled E grid.
  const Temperature T = p.GetCell().temperature;
  const ContinuousEnergy E_s_b_lo =
      F_hi_i != 0 ? EvaluateBeta(E_s_i, F_hi_i - 1, T)
                  : -(E_s * 1e6) / (constants::boltzmann * T);
  const ContinuousEnergy E_s_b_hi =
      F_hi_i != beta_cdfs.size() ? EvaluateBeta(E_s_i, F_hi_i, T) : 20;

  // Evalute interpolated value of beta on the sampled E grid (assuming
  // histogram PDF)
  const Beta b_prime =
      E_s_b_lo + (F - F_lo) / (F_hi - F_lo) * (E_s_b_hi - E_s_b_lo);

  // Scale to preserve thresholds at actual incident energy E. The minimum and
  // maximum values of beta are not included in the dataset since their
  // analytical form is known.
  const Beta b_min = -(E * 1e6) / (constants::boltzmann * T);
  const Beta b_max = 20;
  const Beta E_s_b_min = -(E_s * 1e6) / (constants::boltzmann * T);
  const Beta E_s_b_max = 20;
  return b_min + (b_prime - E_s_b_min) * (b_max - b_min) /
                     (E_s_b_max - E_s_b_min);
}

ThermalScattering::Alpha
ThermalScattering::SampleAlpha(Particle& p, const Beta& b) const noexcept {
}
