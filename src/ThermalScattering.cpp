#include "ThermalScattering.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "Particle.hpp"
#include "ScalarField.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>

// ThermalScattering

//// public

ThermalScattering::ThermalScattering(const pugi::xml_node& tsl_node) noexcept
    : majorant{HDF5DataSet<1>{tsl_node.attribute("majorant").as_string()}
                   .ToContinuousMap()},
      scatter_xs_T{tsl_node.attribute("total_T").as_string()},
      scatter_xs_S{tsl_node.attribute("total_S").as_string()},
      scatter_xs_E{tsl_node.attribute("total_E").as_string()},
      // IIFE
      beta_partitions{[&tsl_node]() noexcept {
        std::vector<BetaPartition> result;
        for (const auto& partition_node : tsl_node.child("beta_partitions")) {
          result.emplace_back(partition_node);
        }
        return result;
      }()},
      // IIFE
      Es{[this]() noexcept {
        std::vector<ContinuousEnergy> result;
        for (const auto& partition : beta_partitions) {
          const auto& Es = partition.E_T_modes.GetAxis(0);
          result.insert(result.end(), Es.cbegin(), Es.cend());
        }
        // incident energies must be monotonically increasing
        ContinuousEnergy prev_E = 0;
        for (const auto& E : result) {
          assert(prev_E < E);
          prev_E = E;
        }
        return result;
      }()},
      // IIFE
      beta_partition_E_ends{[this]() noexcept {
        std::vector<size_t> result;
        size_t current_end = 0;
        for (const auto& partition : beta_partitions) {
          current_end += partition.E_T_modes.GetAxis(0).size();
          result.emplace_back(current_end);
        }
        return result;
      }()},
      // IIFE
      alpha_partitions{[&tsl_node]() noexcept {
        std::vector<AlphaPartition> result;
        for (const auto& partition_node : tsl_node.child("alpha_partitions")) {
          result.emplace_back(partition_node);
        }
        return result;
      }()},
      // IIFE
      betas{[this]() noexcept {
        std::vector<ThermalScattering::Beta> result;
        for (const auto& partition : alpha_partitions) {
          const auto& betas = partition.beta_T_modes.GetAxis(0);
          result.insert(result.end(), betas.cbegin(), betas.cend());
        }
        // betas must be monotonically increasing
        ThermalScattering::Beta prev_beta = 0;
        for (const auto& beta : result) {
          assert(prev_beta < beta);
          prev_beta = beta;
        }
        return result;
      }()},
      // IIFE
      alpha_partition_beta_ends{[this]() {
        std::vector<size_t> result;
        size_t current_end = 0;
        for (const auto& partition : alpha_partitions) {
          current_end += partition.beta_T_modes.GetAxis(0).size();
          result.emplace_back(current_end);
        }
        return result;
      }()},
      beta_cutoff{tsl_node.attribute("beta_cutoff").as_double()},
      alpha_cutoff{tsl_node.attribute("alpha_cutoff").as_double()},
      awr{tsl_node.parent().parent().parent().attribute("awr").as_double()} {}

bool ThermalScattering::IsValid(const Particle& p) const noexcept {
  const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
  return p.GetType() == Particle::Type::neutron && E < cutoff_energy;
}

MicroscopicCrossSection
ThermalScattering::GetMajorant(const Particle& p) const noexcept {
  return majorant.at(std::get<ContinuousEnergy>(p.GetEnergy()));
}

MicroscopicCrossSection
ThermalScattering::GetTotal(const Particle& p) const noexcept {

  // Find index of Energy above and below Particle Energy
  const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
  const auto& Es = scatter_xs_E.GetAxis(0);
  const size_t E_hi_i =
      std::distance(Es.cbegin(), std::upper_bound(Es.cbegin(), Es.cend(), E));
  // We require that E is strictly less than cutoff energy
  assert(E_hi_i != Es.size());
  // If E < E_min, we use the same index for the lower value of E
  const auto below_E_min = E_hi_i == 0;
  const size_t E_lo_i = below_E_min ? E_hi_i : E_hi_i - 1;

  // Find index of Temperature above and below target Temperature
  const Temperature T = p.GetCell().temperature->at(p.GetPosition());
  const auto& Ts = scatter_xs_T.GetAxis(0);
  const size_t candidate_T_hi_i =
      std::distance(Ts.cbegin(), std::upper_bound(Ts.cbegin(), Ts.cend(), T));
  // If T >= T_max, we use T_max for the upper value
  const auto above_T_max = candidate_T_hi_i == Ts.size();
  const size_t T_hi_i = above_T_max ? candidate_T_hi_i - 1 : candidate_T_hi_i;
  // If T < T_min, we use T_min as the lower value
  const auto below_T_min = T_hi_i == 0;
  const size_t T_lo_i = below_T_min ? T_hi_i : T_hi_i - 1;

  // Evaluate neighboring points
  const auto xs_E_lo_T_lo = EvaluateInelastic(E_lo_i, T_lo_i);
  const auto xs_E_lo_T_hi = EvaluateInelastic(E_lo_i, T_hi_i);
  const auto xs_E_hi_T_lo = EvaluateInelastic(E_hi_i, T_lo_i);
  const auto xs_E_hi_T_hi = EvaluateInelastic(E_hi_i, T_hi_i);

  // Linear interpolation in energy
  const auto E_lo = Es.at(E_lo_i);
  const auto E_hi = Es.at(E_hi_i);
  const auto r_E = below_E_min ? 1 : (E - E_lo) / (E_hi - E_lo);
  const auto xs_T_lo = xs_E_lo_T_lo + r_E * (xs_E_hi_T_lo - xs_E_lo_T_lo);
  const auto xs_T_hi = xs_E_lo_T_hi + r_E * (xs_E_hi_T_hi - xs_E_lo_T_hi);

  // Linear interpolation in temperature
  const auto T_lo = Ts.at(T_lo_i);
  const auto T_hi = Ts.at(T_hi_i);
  const auto r_T = below_T_min   ? 1
                   : above_T_max ? 0
                                 : (T - T_lo) / (T_hi - T_lo);
  return xs_T_lo + r_T * (xs_T_hi - xs_T_lo);
}

void ThermalScattering::Scatter(Particle& p) const {
  // get incident energy and target temperature
  const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
  const Temperature T = p.GetCell().temperature->at(p.GetPosition());
  // sample beta then alpha
  const auto beta = SampleBeta(p, E, T);
  const auto alpha = SampleAlpha(p, beta, E, T);
  // convert to outgoing energy and cosine
  const auto E_p = E + beta * constants::boltzmann * T;
  const auto mu = (E + E_p - alpha * awr * constants::boltzmann * T) /
                  (2 * std::sqrt(E * E_p));
  p.Scatter(mu, E_p);
}

//// private

// BetaPartition

ThermalScattering::BetaPartition::BetaPartition(
    const pugi::xml_node& partition_node)
    : CDF_modes{partition_node.attribute("CDF").as_string()},
      singular_values{partition_node.attribute("S").as_string()},
      E_T_modes{partition_node.attribute("E_T").as_string()} {}

ThermalScattering::Beta ThermalScattering::BetaPartition::Evaluate(
    const size_t cdf_index, const size_t E_index,
    Temperature T) const noexcept {

  // Find index of Temperature above and below target Temperature
  const auto& Ts = E_T_modes.GetAxis(1);
  const size_t candidate_T_hi_i =
      std::distance(Ts.cbegin(), std::upper_bound(Ts.cbegin(), Ts.cend(), T));
  // If T >= T_max, we use T_max for the upper value
  const auto above_T_max = candidate_T_hi_i == Ts.size();
  const size_t T_hi_i = above_T_max ? candidate_T_hi_i - 1 : candidate_T_hi_i;
  // If T < T_min, we use T_min as the lower value
  const auto below_T_min = T_hi_i == 0;
  const size_t T_lo_i = below_T_min ? T_hi_i : T_hi_i - 1;

  // Reconstruct betas
  const size_t max_order = singular_values.GetAxis(0).size();

  // handle single temperature case
  Beta beta_T_hi = 0;
  for (size_t order = 0; order < max_order; order++) {
    beta_T_hi += singular_values.at(order) * CDF_modes.at(cdf_index, order) *
                 E_T_modes.at(E_index, T_hi_i, order);
  }
  if (T_hi_i == T_lo_i) {
    return beta_T_hi;
  }

  // the multiple temperature case
  Beta beta_T_lo = 0;
  for (size_t order = 0; order < max_order; order++) {
    beta_T_lo += singular_values.at(order) * CDF_modes.at(cdf_index, order) *
                 E_T_modes.at(E_index, T_lo_i, order);
  }

  // Linearly interpolate
  const Temperature T_hi = Ts.at(T_hi_i);
  const Temperature T_lo = Ts.at(T_lo_i);
  const Real interpolated_beta =
      beta_T_lo + (beta_T_hi - beta_T_lo) / (T_hi - T_lo) * (T - T_lo);
  return interpolated_beta;
}

// AlphaPartition

ThermalScattering::AlphaPartition::AlphaPartition(
    const pugi::xml_node& partition_node)
    : CDF_modes{partition_node.attribute("CDF").as_string()},
      singular_values{partition_node.attribute("S").as_string()},
      beta_T_modes{partition_node.attribute("beta_T").as_string()} {}

ThermalScattering::Alpha ThermalScattering::AlphaPartition::Evaluate(
    const size_t cdf_index, const size_t beta_index, Temperature T) const {

  // Find index of Temperature above and below target Temperature
  const auto& Ts = beta_T_modes.GetAxis(1);
  const size_t candidate_T_hi_i =
      std::distance(Ts.cbegin(), std::upper_bound(Ts.cbegin(), Ts.cend(), T));
  // If T >= T_max, we use T_max for the upper value
  const auto above_T_max = candidate_T_hi_i == Ts.size();
  const size_t T_hi_i = above_T_max ? candidate_T_hi_i - 1 : candidate_T_hi_i;
  // If T < T_min, we use T_min as the lower value
  const auto below_T_min = T_hi_i == 0;
  const size_t T_lo_i = below_T_min ? T_hi_i : T_hi_i - 1;

  // Reconstruct alphas
  const size_t max_order = singular_values.GetAxis(0).size();

  // handle single temperature case
  Alpha alpha_T_hi = 0;
  for (size_t order = 0; order < max_order; order++) {
    alpha_T_hi += singular_values.at(order) * CDF_modes.at(cdf_index, order) *
                  beta_T_modes.at(beta_index, T_hi_i, order);
  }
  if (T_hi_i == T_lo_i) {
    return alpha_T_hi;
  }

  // the multiple temperature case
  Alpha alpha_T_lo = 0;
  for (size_t order = 0; order < max_order; order++) {
    alpha_T_lo += singular_values.at(order) * CDF_modes.at(cdf_index, order) *
                  beta_T_modes.at(beta_index, T_lo_i, order);
  }

  // Linearly interpolate
  const Temperature T_hi = Ts.at(T_hi_i);
  const Temperature T_lo = Ts.at(T_lo_i);
  const Real interpolated_alpha =
      alpha_T_lo + (alpha_T_hi - alpha_T_lo) / (T_hi - T_lo) * (T - T_lo);
  return interpolated_alpha;
}

// ThermalScattering

MacroscopicCrossSection ThermalScattering::EvaluateInelastic(
    const size_t E_index, const size_t T_index) const noexcept {
  const size_t max_order = scatter_xs_S.GetAxis(0).size();
  MicroscopicCrossSection result = 0;
  for (size_t order = 0; order < max_order; order++) {
    result += scatter_xs_S.at(order) * scatter_xs_E.at(E_index, order) *
              scatter_xs_T.at(T_index, order);
  }
  return result;
}

ThermalScattering::Beta ThermalScattering::SampleBeta(
    Particle& p, const ContinuousEnergy E, const Temperature T) const {

  // find index of energy value strictly greater than E
  const size_t E_hi_i =
      std::distance(Es.cbegin(), std::upper_bound(Es.cbegin(), Es.cend(), E));
  // We require that E is strictly less than cutoff energy
  assert(E_hi_i != Es.size());

  // Determine interpolation factor. If E < E_min, we force the use of the grid
  // at E_min by setting r = 1. This way, index of the sampled energy E_s_i is
  // always valid.
  const Real r = E_hi_i != 0 ? (E - Es.at(E_hi_i - 1)) /
                                   (Es.at(E_hi_i) - Es.at(E_hi_i - 1))
                             : 1;
  const size_t E_s_i =
      r <= std::uniform_real_distribution{}(p.rng) ? E_hi_i - 1 : E_hi_i;
  const ContinuousEnergy E_s = Es.at(E_s_i);

  // Get partition which contains the sampled energy and get sampled energy
  // index in that partition
  const size_t P_s_i = std::distance(
      beta_partition_E_ends.cbegin(),
      std::upper_bound(
          beta_partition_E_ends.cbegin(), beta_partition_E_ends.cend(), E_s_i));
  const auto& P_s = beta_partitions.at(P_s_i);
  // TODO: Change beta_partition_E_ends to beta_partition_E_begins to avoid
  // this if statement
  const size_t E_s_i_local =
      P_s_i == 0 ? E_s_i : E_s_i - beta_partition_E_ends.at(P_s_i - 1);
  // sample a CDF value
  const Real F = std::uniform_real_distribution{}(p.rng);
  // find index of CDF value strictly greater than sampled CDF value
  const auto& Fs = P_s.CDF_modes.GetAxis(0);
  const size_t F_hi_i =
      std::distance(Fs.cbegin(), std::upper_bound(Fs.cbegin(), Fs.cend(), F));

  // Evaluate nearest Fs on the sampled E grid.
  const Real F_lo = F_hi_i != 0 ? Fs.at(F_hi_i - 1) : 0;
  const Real F_hi = F_hi_i != Fs.size() ? Fs.at(F_hi_i) : 1;

  // Evaluate nearest betas on the sampled E grid.
  const auto b_s_min = -E_s / (constants::boltzmann * T);
  const ContinuousEnergy E_s_b_lo =
      F_hi_i != 0 ? P_s.Evaluate(F_hi_i - 1, E_s_i_local, T) : b_s_min;
  const ContinuousEnergy E_s_b_hi = F_hi_i != Fs.size()
                                        ? P_s.Evaluate(F_hi_i, E_s_i_local, T)
                                        : beta_cutoff;

  // Evaluate interpolated value of beta on the sampled E grid (assuming
  // histogram PDF)
  const auto b_prime =
      E_s_b_lo + (F - F_lo) / (F_hi - F_lo) * (E_s_b_hi - E_s_b_lo);

  // TODO: fix strange behavior near b_prime = 0 while preserving thresholds
  const auto b_min = -E / (constants::boltzmann * T);
  return b_min +
         (b_prime - b_s_min) / (beta_cutoff - b_s_min) * (beta_cutoff - b_min);
}

ThermalScattering::Alpha ThermalScattering::SampleAlpha(
    Particle& p, const Beta& b, ContinuousEnergy E, Temperature T) const {

  // assume S(a,b) = S(a,-b)
  const Beta abs_b = std::abs(b);
  // get sign of beta: https://stackoverflow.com/a/4609795/5101335
  const int_fast8_t sgn_b = (0 < b) - (b < 0);
  // find index of beta value strictly greater than abs_b
  const size_t b_hi_i = std::distance(
      betas.cbegin(), std::upper_bound(betas.cbegin(), betas.cend(), abs_b));

  // For any given abs_b, we have two beta grids to choose from:
  // b_lo <= abs_b < b_hi
  // For positive b, if beta_cutoff <= b_hi, we force the use of b_lo
  // For negative b, if -b_hi < - E / (k * T), we force the use of b_lo
  const bool snap_to_lower =
      ((sgn_b == 1 && beta_cutoff <= betas.at(b_hi_i)) ||
       (sgn_b == -1 && -betas.at(b_hi_i) < -E / (constants::boltzmann * T)));
  // Handle case where abs_b is less than the smallest beta provided, b_min
  const bool snap_to_min = b_hi_i == 0;
  // Determine interpolation factor and sample a beta grid to use, b_s
  const Real r = snap_to_lower ? 0
                 : snap_to_min
                     ? 1
                     : (abs_b - betas.at(b_hi_i - 1) /
                                    (betas.at(b_hi_i) - betas.at(b_hi_i - 1)));
  const size_t b_s_i =
      r <= std::uniform_real_distribution{}(p.rng) ? b_hi_i - 1 : b_hi_i;
  Beta b_s = sgn_b * betas.at(b_s_i);

  // compute minimum and maximum possible value of alpha
  const auto sqrt_E = std::sqrt(E);
  const auto b_s_sqrt_E_bkT = std::sqrt(E + b_s * constants::boltzmann * T);
  const auto akT = awr * constants::boltzmann * T;
  // minimum and maximum alpha values at sampled beta
  const Alpha b_s_a_min = std::pow(sqrt_E - b_s_sqrt_E_bkT, 2) / akT;
  const Alpha b_s_a_max = std::pow(sqrt_E + b_s_sqrt_E_bkT, 2) / akT;
  assert(b_s_a_max < alpha_cutoff);

  // Get partition which contains the sampled beta and get sampeld beta index
  // in that partition
  const size_t P_s_i = std::distance(
      alpha_partition_beta_ends.cbegin(),
      std::upper_bound(
          alpha_partition_beta_ends.cbegin(), alpha_partition_beta_ends.cend(),
          b_s_i));
  const auto& P_s = alpha_partitions.at(P_s_i);
  // TODO: Change alpha_partition_beta_ends to alpha_partition_beta_begins to
  // avoid this if statement
  const size_t b_s_i_local =
      P_s_i == 0 ? b_s_i : b_s_i - alpha_partition_beta_ends.at(P_s_i - 1);

  const auto& Fs = P_s.CDF_modes.GetAxis(0);
  // helper function to find the CDF value that would return `a`. (assuming
  // histogram PDF)
  const auto find_cdf = [this, &Fs, &P_s, &b_s_i_local,
                         &T](const Alpha& a) -> Real {
    // find index of CDF value that evaluates to be strictly greater than `a`
    const size_t F_a_hi_i = std::distance(
        Fs.cbegin(),
        std::upper_bound(
            Fs.cbegin(), Fs.cend(), a,
            [&Fs, &P_s, &b_s_i_local, &T](const Alpha& a, const Real& cdf) {
              // https://en.cppreference.com/w/cpp/language/operator_arithmetic
              const size_t cdf_i = &cdf - &Fs.front();
              return a < P_s.Evaluate(cdf_i, b_s_i_local, T);
            }));
    // evaluate nearest CDFs
    const Real F_a_lo = F_a_hi_i != 0 ? Fs.at(F_a_hi_i - 1) : 0;
    const Real F_a_hi = F_a_hi_i != Fs.size() ? Fs.at(F_a_hi_i) : 1;
    // evaluate nearest alphas
    const Alpha a_lo =
        F_a_hi_i != 0 ? P_s.Evaluate(F_a_hi_i - 1, b_s_i_local, T) : 0;
    const Alpha a_hi = F_a_hi_i != Fs.size()
                           ? P_s.Evaluate(F_a_hi_i, b_s_i_local, T)
                           : alpha_cutoff;
    // assume histogram PDF
    return F_a_lo + (a - a_lo) * (F_a_hi - F_a_lo) / (a_hi - a_lo);
  };
  // CDFs corresponding to b_s_a_min and b_s_a_max
  const auto F_min = find_cdf(b_s_a_min);
  const auto F_max = find_cdf(b_s_a_max);

  // Keep sampling until a physical value of alpha is obtained
  for (size_t resamples = 0; resamples < constants::alpha_resample_limit;
       resamples++) {
    // sample a CDF value which is scaled to return a result in [b_s_a_min,
    // b_s_a_max]
    const Real F =
        F_min + std::uniform_real_distribution{}(p.rng) * (F_max - F_min);
    // find index of CDF value strictly greater than sampled CDF value
    const size_t F_hi_i =
        std::distance(Fs.cbegin(), std::upper_bound(Fs.cbegin(), Fs.cend(), F));

    // Evaluate nearest CDFs on the sampled beta grid
    const Real F_lo = F_hi_i != 0 ? Fs.at(F_hi_i - 1) : 0;
    const Real F_hi = F_hi_i != Fs.size() ? Fs.at(F_hi_i) : 1;

    // Evaluate nearest alphas on the sampled beta grid.
    const Real b_s_a_lo =
        F_hi_i != 0 ? P_s.Evaluate(F_hi_i - 1, b_s_i_local, T) : 0;
    const Real b_s_a_hi = F_hi_i != Fs.size()
                              ? P_s.Evaluate(F_hi_i, b_s_i_local, T)
                              : alpha_cutoff;

    // Evalute interpolated value of alpha on the sampled beta grid (assuming
    // histogram PDF)
    const auto a_prime =
        b_s_a_lo + (F - F_lo) / (F_hi - F_lo) * (b_s_a_hi - b_s_a_lo);

    if (b_s_a_min < a_prime && a_prime < b_s_a_max) {
      // Scale to preserve thresholds at actual beta value b. The minimum and
      // maximum values of alpha have analytical forms which are known.
      const auto b_sqrt_E_bkT = std::sqrt(E + b * constants::boltzmann * T);
      const Alpha b_a_min = std::pow(sqrt_E - b_sqrt_E_bkT, 2) / akT;
      const Alpha b_a_max = std::pow(sqrt_E + b_sqrt_E_bkT, 2) / akT;
      return b_a_min + (a_prime - b_s_a_min) * (b_a_max - b_a_min) /
                           (b_s_a_max - b_s_a_min);
    }
  }
  throw std::runtime_error(
      "Number of thermal scattering alpha resamples exceeded limit: " +
      std::to_string(constants::alpha_resample_limit));
}
