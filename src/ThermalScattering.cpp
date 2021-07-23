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
  // Data is only provided at discrete energy and CDF values.
  ContinuousEnergy lo_E; // nearest lower value of energy found
  ContinuousEnergy hi_E; // nearest higher value of energy found
  Real lo_F;             // nearest lower value of CDF found
  Real hi_F;             // nearest higher value of CDf found
  // four nearest values of beta are needed for bilinear interpolation
  Beta lo_E_lo_F; // beta at lower value of energy, lower value of CDF
  Beta lo_E_hi_F; // beta at lower value of energy, higher value of CDF
  Beta hi_E_lo_F; // beta at higher value of energy, lower value of CDF
  Beta hi_E_hi_F; // beta at higher value of energy, higher value of CDF

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
      // Find nearest E and CDF values
      hi_E = beta_energies.at(hi_E_i);
      lo_E = beta_energies.at(hi_E_i - 1);
      lo_F = beta_cdfs.at(hi_F_i - 1);
      hi_F = beta_cdfs.at(hi_F_i);
      // Evaluate beta at nearest points
      hi_E_lo_F = EvaluateBeta(hi_E_i, hi_F_i - 1, T);
      hi_E_hi_F = EvaluateBeta(hi_E_i, hi_F_i, T);
      lo_E_lo_F = EvaluateBeta(hi_E_i - 1, hi_F_i - 1, T);
      lo_E_hi_F = EvaluateBeta(hi_E_i - 1, hi_F_i, T);
    }
    // Less common case, E < E_min
    else {
      // Find nearest E and CDF values
      hi_E = beta_energies.at(hi_E_i);
      lo_E = 0;
      lo_F = beta_cdfs.at(hi_F_i - 1);
      hi_F = beta_cdfs.at(hi_F_i);
      // Evaluate beta at nearest points
      hi_E_lo_F = EvaluateBeta(hi_E_i, hi_F_i - 1, T);
      hi_E_hi_F = EvaluateBeta(hi_E_i, hi_F_i, T);
      lo_E_lo_F = hi_E_lo_F;
      lo_E_hi_F = hi_E_hi_F;
    }
  }
  // Less common case, cdf < cdf_min
  else if (hi_F_i == 0) {
    // More common case, E_min <= E < E_max
    if (hi_E_i != 0) {
      // Find nearest E and CDF values
      hi_E = beta_energies.at(hi_E_i);
      lo_E = beta_energies.at(hi_E_i - 1);
      lo_F = 0;
      hi_F = beta_cdfs.at(hi_F_i);
      // Evaluate beta at nearest points
      hi_E_lo_F = -(beta_energies.at(hi_E_i) * 1e6) / (constants::boltzmann * T);
      hi_E_hi_F = EvaluateBeta(hi_E_i, hi_F_i, T);
      lo_E_lo_F =
          -(beta_energies.at(hi_E_i - 1) * 1e6) / (constants::boltzmann * T);
      lo_E_hi_F = EvaluateBeta(hi_E_i - 1, hi_F_i, T);
    }
    // Less common case, E < E_min
    else {
      // Find nearest E and CDF values
      hi_E = beta_energies.at(hi_E_i);
      lo_E = 0;
      lo_F = 0;
      hi_F = beta_cdfs.at(hi_F_i);
      // Evaluate beta at nearest points
      hi_E_lo_F = -(beta_energies.at(hi_E_i) * 1e6) / (constants::boltzmann * T);
      hi_E_hi_F = EvaluateBeta(hi_E_i, hi_F_i, T);
      lo_E_lo_F = 0;         // - E / (k * T) = 0
      lo_E_hi_F = hi_E_hi_F; // use same value as hi_E
    }
  }
  // Least common case, cdf_max <= cdf
  // implicit condition: hi_F_i == beta_cdfs.size()
  else {
    // More common case, E_min <= E < E_max
    if (hi_E_i != 0) {
      // Find nearest E and CDF values
      hi_E = beta_energies.at(hi_E_i);
      lo_E = beta_energies.at(hi_E_i - 1);
      lo_F = beta_cdfs.at(hi_F_i - 1);
      hi_F = 1;
      // Evaluate beta at nearest points
      hi_E_lo_F = EvaluateBeta(hi_E_i, hi_F_i - 1, T);
      hi_E_hi_F = 20;
      lo_E_lo_F = EvaluateBeta(hi_E_i - 1, hi_F_i - 1, T);
      lo_E_hi_F = 20;
    }
    // Less common case, E < E_min
    else {
      // Find nearest E and CDF values
      hi_E = beta_energies.at(hi_E_i);
      lo_E = 0;
      lo_F = beta_cdfs.at(hi_F_i - 1);
      hi_F = 1;
      // Evaluate beta at nearest points
      hi_E_lo_F = EvaluateBeta(hi_E_i, hi_F_i - 1, T);
      hi_E_hi_F = 20;
      lo_E_lo_F = hi_E_lo_F; // use same value as hi_E
      lo_E_hi_F = 20;
    }
  }
  // bilinear interpolation in E and CDF
  return (lo_E_lo_F * (hi_E - E) * (hi_F - F) +
          hi_E_lo_F * (E - lo_E) * (hi_F - F) +
          lo_E_hi_F * (hi_E - E) * (F - lo_F) +
          hi_E_hi_F * (E - lo_E) * (F - lo_F)) /
         ((hi_E - lo_E) * (hi_F - lo_F));
}

ThermalScattering::Alpha
ThermalScattering::SampleAlpha(Particle& p, const Beta& b) const noexcept {
  // Data is only provided at discrete beta and CDF values.
  Beta lo_b; // nearest lower value of energy found
  Beta hi_b; // nearest higher value of energy found
  Real lo_F; // nearest lower value of CDF found
  Real hi_F; // nearest higher value of CDf found
  // four nearest values of alpha are needed for bilinear interpolation
  Alpha lo_b_lo_F; // alpha at lower value of beta, lower value of CDF
  Alpha lo_b_hi_F; // alpha at lower value of beta, higher value of CDF
  Alpha hi_b_lo_F; // alpha at higher value of beta, lower value of CDF
  Alpha hi_b_hi_F; // alpha at higher value of beta, higher value of CDF

  // temperature of Cell occupied by Particle
  const Temperature T = p.GetCell().temperature;
  // find index of beta value strictly greater than b
  const size_t hi_b_i = std::distance(
      alpha_betas.cbegin(),
      std::upper_bound(alpha_betas.cbegin(), alpha_betas.cend(), b));
  // sample a CDF value
  const Real F = std::uniform_real_distribution{}(p.rng);
  // find index of CDF value strictly greater than sampled CDF value
  const size_t hi_F_i = std::distance(
      beta_cdfs.cbegin(),
      std::upper_bound(beta_cdfs.cbegin(), beta_cdfs.cend(), F));

  // In all the following cases, we require that a datapoint for b_min = 0
  // exists
  assert(hi_b_i != 0);
  // In the following comments, cdf_min and cdf_max are the smallest and
  // greatest CDF values that appear in the dataset, respectively.
  // b_max is the greatest beta value that appears in the dataset.

  // The most common case, cdf_min <= cdf < cdf_max
  if (hi_F_i != 0 && hi_F_i != alpha_cdfs.size()) {
    // The most common case, 0 <= b < b_max
    if (hi_b_i != alpha_betas.size()) {
      // Find nearest beta and CDF values
      lo_b = alpha_betas.at(hi_b_i - 1);
      hi_b = alpha_betas.at(hi_b_i);
      lo_F = alpha_cdfs.at(hi_F_i - 1);
      hi_F = alpha_cdfs.at(hi_F_i);
      // Evaluate alpha at nearest points
      lo_b_lo_F = EvaluateAlpha(hi_b_i - 1, hi_F_i - 1, T);
      lo_b_hi_F = EvaluateAlpha(hi_b_i - 1, hi_F_i, T);
      hi_b_lo_F = EvaluateAlpha(hi_b_i, hi_F_i - 1, T);
      hi_b_hi_F = EvaluateAlpha(hi_b_i, hi_F_i, T);
    }
    // Less common case, b_max <= b
    else {
      // Find nearest beta and CDF values
      lo_b = alpha_betas.at(hi_b_i - 1);
      hi_b = lo_b;
      lo_F = alpha_cdfs.at(hi_F_i - 1);
      hi_F = alpha_cdfs.at(hi_F_i);
      // Evaluate alpha at nearest points
      lo_b_lo_F = EvaluateAlpha(hi_b_i - 1, hi_F_i - 1, T);
      lo_b_hi_F = EvaluateAlpha(hi_b_i - 1, hi_F_i, T);
      hi_b_lo_F = lo_b_lo_F;
      hi_b_hi_F = lo_b_hi_F;
    }
  }
  // Less common case, cdf < cdf_min
  else if (hi_F_i == 0) {
    // More common case, 0 <= b < b_max
    if (hi_b_i != alpha_betas.size()) {
      // Find nearest beta and CDF values
      lo_b = alpha_betas.at(hi_b_i - 1);
      hi_b = alpha_betas.at(hi_b_i);
      lo_F = 0;
      hi_F = alpha_cdfs.at(hi_F_i);
      // Evaluate alpha at nearest points
      const auto E = std::get<ContinuousEnergy>(p.GetEnergy()) * 1e6;
      lo_b_lo_F =
          std::pow(
              std::sqrt(E) - std::sqrt(E + lo_b * constants::boltzmann * T),
              2) /
          (awr * constants::boltzmann * T);
      lo_b_hi_F = EvaluateBeta(hi_b_i - 1, hi_F_i, T);
      hi_b_lo_F =
          std::pow(
              std::sqrt(E) - std::sqrt(E + hi_b * constants::boltzmann * T),
              2) /
          (awr * constants::boltzmann * T);
      hi_b_hi_F = EvaluateBeta(hi_b_i, hi_F_i, T);
    }
    // Less common case, b_max <= b
    else {
      // Find nearest beta and CDF values
      lo_b = alpha_betas.at(hi_b_i - 1);
      hi_b = lo_b;
      lo_F = 0;
      hi_F = alpha_cdfs.at(hi_F_i);
      // Evaluate alpha at nearest points
      const auto E = std::get<ContinuousEnergy>(p.GetEnergy()) * 1e6;
      lo_b_lo_F =
          std::pow(
              std::sqrt(E) - std::sqrt(E + lo_b * constants::boltzmann * T),
              2) /
          (awr * constants::boltzmann * T);
      lo_b_hi_F = EvaluateAlpha(hi_b_i - 1, hi_F_i, T);
      hi_b_lo_F = lo_b_lo_F;
      hi_b_hi_F = lo_b_hi_F;
    }
  }
  // Less common case, cdf_max <= cdf
  // implicit condition: hi_F_i == beta_cdfs.size()
  else {
    // More common case, E_min <= b < E_max
    if (hi_b_i != alpha_betas.size()) {
      // Find nearest beta and CDF values
      lo_b = alpha_betas.at(hi_b_i - 1);
      hi_b = alpha_betas.at(hi_b_i);
      lo_F = alpha_cdfs.at(hi_F_i - 1);
      hi_F = 1;
      // Evaluate alpha at nearest points
      const auto E = std::get<ContinuousEnergy>(p.GetEnergy()) * 1e6;
      lo_b_lo_F = EvaluateAlpha(hi_b_i - 1, hi_F_i - 1, T);
      lo_b_hi_F =
          std::pow(
              std::sqrt(E) + std::sqrt(E + lo_b * constants::boltzmann * T),
              2) /
          (awr * constants::boltzmann * T);
      hi_b_lo_F = EvaluateAlpha(hi_b_i, hi_F_i - 1, T);
      hi_b_hi_F =
          std::pow(
              std::sqrt(E) + std::sqrt(E + hi_b * constants::boltzmann * T),
              2) /
          (awr * constants::boltzmann * T);
    }
    // Less common case, b_max <= b
    else {
      // Find nearest beta and CDF values
      lo_b = alpha_betas.at(hi_b_i - 1);
      hi_b = lo_b;
      lo_F = beta_cdfs.at(hi_F_i - 1);
      hi_F = 1;
      // Evaluate alpha at nearest points
      const auto E = std::get<ContinuousEnergy>(p.GetEnergy()) * 1e6;
      lo_b_lo_F = EvaluateBeta(hi_b_i - 1, hi_F_i -1, T);
      lo_b_hi_F =
          std::pow(
              std::sqrt(E) + std::sqrt(E + lo_b * constants::boltzmann * T),
              2) /
          (awr * constants::boltzmann * T);
      hi_b_lo_F = lo_b_lo_F;
      hi_b_hi_F = lo_b_hi_F;
    }
  }
  // bilinear interpolation in beta and CDF
  return (lo_b_lo_F * (hi_b - b) * (hi_F - F) +
          hi_b_lo_F * (b - lo_b) * (hi_F - F) +
          lo_b_hi_F * (hi_b - b) * (F - lo_F) +
          hi_b_hi_F * (b - lo_b) * (F - lo_F)) /
         ((hi_b - lo_b) * (hi_F - lo_F));
}
