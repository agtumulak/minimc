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
ThermalScattering::SampleBetaHistogramPDF(Particle& p) const noexcept {

  // temperature of Cell occupied by Particle
  const Temperature T = p.GetCell().temperature;

  // incident energy of neutron
  const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
  // find index of energy value strictly greater than E
  const size_t hi_E_i = std::distance(
      beta_energies.cbegin(),
      std::upper_bound(beta_energies.cbegin(), beta_energies.cend(), E));

  // sample a CDF value
  const Real F = std::uniform_real_distribution{}(p.rng);
  // find index of CDF value strictly greater than sampled CDF value
  const size_t hi_F_i = std::distance(
      beta_cdfs.cbegin(),
      std::upper_bound(beta_cdfs.cbegin(), beta_cdfs.cend(), F));

  // In all the following cases, we require that E is strictly less than
  // cutoff energy
  assert (hi_E_i != beta_energies.size());
  // In the following comments, cdf_min and cdf_max are the smallest and
  // greatest CDF values that appear in the dataset, respectively.
  // E_min and E_max are the smallest and greatest energy values that appear
  // in the dataset, respectively.

  // The most common case, cdf_min <= cdf < cdf_max
  if (hi_F_i != 0 && hi_F_i != beta_cdfs.size()) {
    // The most common case, E_min <= E < E_max
    if (hi_E_i != 0) {
      const ContinuousEnergy hi_E = beta_energies.at(hi_E_i);
      const ContinuousEnergy lo_E = beta_energies.at(hi_E_i - 1);
      const Real r = (E - lo_E) / (hi_E - lo_E);
      const size_t sampled_E_i =
          r < std::uniform_real_distribution{}(p.rng) ? hi_E_i - 1 : hi_E_i;
      const ContinuousEnergy sampled_E = beta_energies.at(sampled_E_i);
      const Real hi_F = beta_cdfs.at(hi_F_i);
      const Real lo_F = beta_cdfs.at(hi_F_i - 1);
      const Beta sampled_E_hi_b = EvaluateBeta(sampled_E_i, hi_F_i, T);
      const Beta sampled_E_lo_b = EvaluateBeta(sampled_E_i, hi_F_i - 1, T);
      // Sample beta from sampled energy grid using histogram PDF
      // interpolation. Unscaled.
      const Beta b_prime =
          sampled_E_lo_b +
          (F - lo_F) / (hi_F - lo_F) * (sampled_E_hi_b - sampled_E_lo_b);
      // Scale to preserve thresholds. The lower and upper values of beta are
      // not included in the data set since their analytical form is known.
      const Beta b_min = - (E * 1e6) / (constants::boltzmann * T);
      const Beta b_max = 20;
      const Beta sampled_E_b_min =
          -(sampled_E * 1e6) / (constants::boltzmann * T);
      const Beta sampled_E_b_max = 20;
      return b_min + (b_prime - sampled_E_b_min) * (b_max - b_min) /
                         (sampled_E_b_max - sampled_E_b_min);
    }
    // Less common case, E < E_min
    else {
      const ContinuousEnergy hi_E = beta_energies.at(hi_E_i);
      const ContinuousEnergy lo_E = 0; // XXX: Different
      const Real r = 1; // XXX: Different
      const size_t sampled_E_i =
          r < std::uniform_real_distribution{}(p.rng) ? hi_E_i - 1 : hi_E_i;
      const ContinuousEnergy sampled_E = beta_energies.at(sampled_E_i);
      const Real hi_F = beta_cdfs.at(hi_F_i);
      const Real lo_F = beta_cdfs.at(hi_F_i - 1);
      const Beta sampled_E_hi_b = EvaluateBeta(sampled_E_i, hi_F_i, T);
      const Beta sampled_E_lo_b = EvaluateBeta(sampled_E_i, hi_F_i - 1, T);
      // Sample beta from sampled energy grid using histogram PDF
      // interpolation. Unscaled.
      const Beta b_prime =
          sampled_E_lo_b +
          (F - lo_F) / (hi_F - lo_F) * (sampled_E_hi_b - sampled_E_lo_b);
      // Scale to preserve thresholds. The lower and upper values of beta are
      // not included in the data set since their analytical form is known.
      const Beta b_min = - (E * 1e6) / (constants::boltzmann * T);
      const Beta b_max = 20;
      const Beta sampled_E_b_min =
          -(sampled_E * 1e6) / (constants::boltzmann * T);
      const Beta sampled_E_b_max = 20;
      return b_min + (b_prime - sampled_E_b_min) * (b_max - b_min) /
                         (sampled_E_b_max - sampled_E_b_min);
    }
  }
  // Less common case, cdf < cdf_min
  else if (hi_F_i == 0) {
    // More common case, E_min <= E < E_max
    if (hi_E_i != 0) {
      const ContinuousEnergy hi_E = beta_energies.at(hi_E_i);
      const ContinuousEnergy lo_E = beta_energies.at(hi_E_i - 1);
      const Real r = (E - lo_E) / (hi_E - lo_E);
      const size_t sampled_E_i =
          r < std::uniform_real_distribution{}(p.rng) ? hi_E_i - 1 : hi_E_i;
      const ContinuousEnergy sampled_E = beta_energies.at(sampled_E_i);
      const Real hi_F = beta_cdfs.at(hi_F_i);
      const Real lo_F = 0; // XXX: Different
      const Beta sampled_E_hi_b = EvaluateBeta(sampled_E_i, hi_F_i, T);
      const Beta sampled_E_lo_b =
          -(sampled_E * 1e6) / (constants::boltzmann * T); // XXX: Different
      // Sample beta from sampled energy grid using histogram PDF
      // interpolation. Unscaled.
      const Beta b_prime =
          sampled_E_lo_b +
          (F - lo_F) / (hi_F - lo_F) * (sampled_E_hi_b - sampled_E_lo_b);
      // Scale to preserve thresholds. The lower and upper values of beta are
      // not included in the data set since their analytical form is known.
      const Beta b_min = - (E * 1e6) / (constants::boltzmann * T);
      const Beta b_max = 20;
      const Beta sampled_E_b_min =
          -(sampled_E * 1e6) / (constants::boltzmann * T);
      const Beta sampled_E_b_max = 20;
      return b_min + (b_prime - sampled_E_b_min) * (b_max - b_min) /
                         (sampled_E_b_max - sampled_E_b_min);
    }
    // Less common case, E < E_min
    else {
      const ContinuousEnergy hi_E = beta_energies.at(hi_E_i);
      const ContinuousEnergy lo_E = 0; // XXX: Different
      const Real r = 1; // XXX: Different
      const size_t sampled_E_i =
          r < std::uniform_real_distribution{}(p.rng) ? hi_E_i - 1 : hi_E_i;
      const ContinuousEnergy sampled_E = beta_energies.at(sampled_E_i);
      const Real hi_F = beta_cdfs.at(hi_F_i);
      const Real lo_F = 0; // XXX: Different
      const Beta sampled_E_hi_b = EvaluateBeta(sampled_E_i, hi_F_i, T);
      const Beta sampled_E_lo_b =
          -(sampled_E * 1e6) / (constants::boltzmann * T); // XXX: Different
      // Sample beta from sampled energy grid using histogram PDF
      // interpolation. Unscaled.
      const Beta b_prime =
          sampled_E_lo_b +
          (F - lo_F) / (hi_F - lo_F) * (sampled_E_hi_b - sampled_E_lo_b);
      // Scale to preserve thresholds. The lower and upper values of beta are
      // not included in the data set since their analytical form is known.
      const Beta b_min = - (E * 1e6) / (constants::boltzmann * T);
      const Beta b_max = 20;
      const Beta sampled_E_b_min =
          -(sampled_E * 1e6) / (constants::boltzmann * T);
      const Beta sampled_E_b_max = 20;
      return b_min + (b_prime - sampled_E_b_min) * (b_max - b_min) /
                         (sampled_E_b_max - sampled_E_b_min);
    }
  }
  // Least common case, cdf_max <= cdf
  // implicit condition: hi_F_i == beta_cdfs.size()
  else {
    // More common case, E_min <= E < E_max
    if (hi_E_i != 0) {
      const ContinuousEnergy hi_E = beta_energies.at(hi_E_i);
      const ContinuousEnergy lo_E = beta_energies.at(hi_E_i - 1);
      const Real r = (E - lo_E) / (hi_E - lo_E);
      const size_t sampled_E_i =
          r < std::uniform_real_distribution{}(p.rng) ? hi_E_i - 1 : hi_E_i;
      const ContinuousEnergy sampled_E = beta_energies.at(sampled_E_i);
      const Real hi_F = 1; // XXX: Different
      const Real lo_F = beta_cdfs.at(hi_F_i - 1);
      const Beta sampled_E_hi_b = 20; // XXX: Different
      const Beta sampled_E_lo_b = EvaluateBeta(sampled_E_i, hi_F_i - 1, T);
      // Sample beta from sampled energy grid using histogram PDF
      // interpolation. Unscaled.
      const Beta b_prime =
          sampled_E_lo_b +
          (F - lo_F) / (hi_F - lo_F) * (sampled_E_hi_b - sampled_E_lo_b);
      // Scale to preserve thresholds. The lower and upper values of beta are
      // not included in the data set since their analytical form is known.
      const Beta b_min = - (E * 1e6) / (constants::boltzmann * T);
      const Beta b_max = 20;
      const Beta sampled_E_b_min =
          -(sampled_E * 1e6) / (constants::boltzmann * T);
      const Beta sampled_E_b_max = 20;
      return b_min + (b_prime - sampled_E_b_min) * (b_max - b_min) /
                         (sampled_E_b_max - sampled_E_b_min);
    }
    // Less common case, E < E_min
    else {
      const ContinuousEnergy hi_E = beta_energies.at(hi_E_i);
      const ContinuousEnergy lo_E = 0; // XXX: Different
      const Real r = 1; // XXX: Different
      const size_t sampled_E_i =
          r < std::uniform_real_distribution{}(p.rng) ? hi_E_i - 1 : hi_E_i;
      const ContinuousEnergy sampled_E = beta_energies.at(sampled_E_i);
      const Real hi_F = beta_cdfs.at(hi_F_i);
      const Real lo_F = beta_cdfs.at(hi_F_i - 1);
      const Beta sampled_E_hi_b = EvaluateBeta(sampled_E_i, hi_F_i, T);
      const Beta sampled_E_lo_b = EvaluateBeta(sampled_E_i, hi_F_i - 1, T);
      // Sample beta from sampled energy grid using histogram PDF
      // interpolation. Unscaled.
      const Beta b_prime =
          sampled_E_lo_b +
          (F - lo_F) / (hi_F - lo_F) * (sampled_E_hi_b - sampled_E_lo_b);
      // Scale to preserve thresholds. The lower and upper values of beta are
      // not included in the data set since their analytical form is known.
      const Beta b_min = - (E * 1e6) / (constants::boltzmann * T);
      const Beta b_max = 20;
      const Beta sampled_E_b_min =
          -(sampled_E * 1e6) / (constants::boltzmann * T);
      const Beta sampled_E_b_max = 20;
      return b_min + (b_prime - sampled_E_b_min) * (b_max - b_min) /
                         (sampled_E_b_max - sampled_E_b_min);
    }
  }
}

ThermalScattering::Alpha
ThermalScattering::SampleAlpha(Particle& p, const Beta& b) const noexcept {}
