#include "ThermalScattering.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "Nuclide.hpp"
#include "Particle.hpp"
#include "Perturbation/IndirectEffect/IndirectEffect.hpp"
#include "Perturbation/IndirectEffect/Visitor.hpp"
#include "ScalarField.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <memory>
#include <random>
#include <type_traits>
#include <variant>

// ThermalScattering

//// public

Real ThermalScattering::BilinearForm(
    const Real x0, const Real x1, const Real a00, const Real a01,
    const Real a10, const Real a11, const Real y0, const Real y1) noexcept {
  return x0 * a00 * y0 + x0 * a01 * y1 + x1 * a10 * y0 + x1 * a11 * y1;
}

ThermalScattering::ThermalScattering(
    const pugi::xml_node& tnsl_node, const Nuclide& target) noexcept
    : majorant{HDF5DataSet<1>{tnsl_node.attribute("majorant").as_string()}
                   .ToContinuousMap()},
      scatter_xs_T{tnsl_node.attribute("total_T").as_string()},
      scatter_xs_S{tnsl_node.attribute("total_S").as_string()},
      scatter_xs_E{tnsl_node.attribute("total_E").as_string()},
      // IIFE
      beta_partitions{[&tnsl_node]() noexcept {
        std::vector<BetaPartition> result;
        for (const auto& partition_node : tnsl_node.child("beta_partitions")) {
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
      alpha_partitions{[&tnsl_node]() noexcept {
        std::vector<AlphaPartition> result;
        for (const auto& partition_node : tnsl_node.child("alpha_partitions")) {
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
      alpha_partition_beta_begins{[this]() {
        std::vector<size_t> result{0};
        for (const auto& partition : alpha_partitions) {
          result.emplace_back(
              result.back() + partition.beta_T_modes.GetAxis(0).size());
        }
        assert(result.back() == betas.size());
        return result;
      }()},
      beta_cutoff{tnsl_node.attribute("beta_cutoff").as_double()},
      alpha_cutoff{tnsl_node.attribute("alpha_cutoff").as_double()},
      target{target} {}

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

void ThermalScattering::Scatter(Particle& p) const noexcept {
  // get incident energy and target temperature
  const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
  const Temperature T = p.GetCell().temperature->at(p.GetPosition());
  // sample beta then alpha
  const auto beta = SampleBeta(p, E, T);
  const auto alpha = SampleAlpha(p, beta, E, T);
  // convert to outgoing energy and cosine
  const auto E_p = E + beta * constants::boltzmann * T;
  const auto mu = (E + E_p - alpha * target.awr * constants::boltzmann * T) /
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
      beta_T_modes{partition_node.attribute("beta_T").as_string()} {
  // The following assumptions about CDF values will be used in sampling
  const auto& Fs = CDF_modes.GetAxis(0);
  assert(Fs.size() > 0);
}

ThermalScattering::Alpha ThermalScattering::AlphaPartition::Evaluate(
    const size_t cdf_index, const size_t local_beta_index,
    const size_t T_index) const {
  Alpha alpha = 0;
  for (size_t order = 0; order < singular_values.GetAxis(0).size(); order++) {
    alpha += singular_values.at(order) * CDF_modes.at(cdf_index, order) *
             beta_T_modes.at(local_beta_index, T_index, order);
  }
  return alpha;
}

Real ThermalScattering::AlphaPartition::FindCDF(
    const Alpha a, const size_t b_s_i_local, const Temperature T,
    const Alpha alpha_cutoff) const noexcept {
  // Evaluate nearest temperatures on the sampled beta grid
  const auto& Ts = beta_T_modes.GetAxis(1);
  // find neighboring temperature value indices
  const size_t T_hi_i =
      std::distance(Ts.cbegin(), std::upper_bound(Ts.cbegin(), Ts.cend(), T));
  // identify two possible edge cases for temperature
  const auto T_hi_exists = T_hi_i < Ts.size();
  const auto T_lo_exists = T_hi_i != 0;
  // find index of CDF value that evaluates to be strictly greater than `a`
  const auto& Fs = CDF_modes.GetAxis(0);
  // TODO: Factor out common terms
  if (T_hi_exists && T_lo_exists) {
    // the most common case
    const size_t F_hi_i = std::distance(
        Fs.cbegin(),
        std::upper_bound(
            Fs.cbegin(), Fs.cend(), a,
            [this, &Fs, &b_s_i_local, &T_hi_i, &Ts,
             &T](const Alpha& a, const Real& F) {
              // https://en.cppreference.com/w/cpp/language/operator_arithmetic
              const size_t F_i = &F - &Fs.front();
              const auto a_T_hi = Evaluate(F_i, b_s_i_local, T_hi_i);
              const auto a_T_lo = Evaluate(F_i, b_s_i_local, T_hi_i - 1);
              // linearly interpolate in temperature
              const auto T_hi = Ts.at(T_hi_i);
              const auto T_lo = Ts.at(T_hi_i - 1);
              return a <
                     a_T_lo + (T - T_lo) / (T_hi - T_lo) * (a_T_hi - a_T_lo);
            }));
    // identify two possible edge cases for CDF
    const auto F_hi_exists = F_hi_i < Fs.size();
    const auto F_lo_exists = F_hi_i != 0;
    if (F_hi_exists && F_lo_exists) {
      // the most common case (1)
      const auto a_F_hi_T_hi = Evaluate(F_hi_i, b_s_i_local, T_hi_i);
      const auto a_F_hi_T_lo = Evaluate(F_hi_i, b_s_i_local, T_hi_i - 1);
      const auto a_F_lo_T_hi = Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i);
      const auto a_F_lo_T_lo = Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i - 1);
      // linearly interpolate upper and lower value of F in temperature
      const auto T_hi = Ts.at(T_hi_i);
      const auto T_lo = Ts.at(T_hi_i - 1);
      const auto r_T = (T - T_lo) / (T_hi - T_lo);
      const auto a_F_hi = a_F_hi_T_lo + r_T * (a_F_hi_T_hi - a_F_hi_T_lo);
      const auto a_F_lo = a_F_lo_T_lo + r_T * (a_F_lo_T_hi - a_F_lo_T_lo);
      // return value of F that would return `a`
      const auto r_a = (a - a_F_lo) / (a_F_hi - a_F_lo);
      const auto F_hi = Fs.at(F_hi_i);
      const auto F_lo = Fs.at(F_hi_i - 1);
      return F_lo + r_a * (F_hi - F_lo);
    }
    else if (F_hi_exists && !F_lo_exists) {
      // use zero alpha at zero CDF (4)
      const auto a_F_hi_T_hi = Evaluate(F_hi_i, b_s_i_local, T_hi_i);
      const auto a_F_hi_T_lo = Evaluate(F_hi_i, b_s_i_local, T_hi_i - 1);
      // linearly interpolate upper value of F in temperature
      const auto T_hi = Ts.at(T_hi_i);
      const auto T_lo = Ts.at(T_hi_i - 1);
      const auto r_T = (T - T_lo) / (T_hi - T_lo);
      const auto a_F_hi = a_F_hi_T_lo + r_T * (a_F_hi_T_hi - a_F_hi_T_lo);
      const Alpha a_F_lo = 0;
      // return value of F that would return `a`
      const auto r_a = (a - a_F_lo) / (a_F_hi - a_F_lo);
      const auto F_hi = Fs.at(F_hi_i);
      const Real F_lo = 0;
      return F_lo + r_a * (F_hi - F_lo);
    }
    else if (!F_hi_exists && F_lo_exists) {
      // use alpha cutoff at unity CDF (5)
      const auto a_F_lo_T_hi = Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i);
      const auto a_F_lo_T_lo = Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i - 1);
      // linearly interpolate lower value of F in temperature
      const auto T_hi = Ts.at(T_hi_i);
      const auto T_lo = Ts.at(T_hi_i - 1);
      const auto r_T = (T - T_lo) / (T_hi - T_lo);
      const auto a_F_hi = alpha_cutoff;
      const auto a_F_lo = a_F_lo_T_lo + r_T * (a_F_lo_T_hi - a_F_lo_T_lo);
      // return value of F that would return `a`
      const auto r_a = (a - a_F_lo) / (a_F_hi - a_F_lo);
      const Real F_hi = 1;
      const auto F_lo = Fs.at(F_hi_i - 1);
      return F_lo + r_a * (F_hi - F_lo);
    }
    else {
      assert(false);
    }
  }
  else if (T_hi_exists && !T_lo_exists) {
    // snap to upper temperature (2)
    const size_t F_hi_i = std::distance(
        Fs.cbegin(),
        std::upper_bound(
            Fs.cbegin(), Fs.cend(), a,
            [this, &Fs, &b_s_i_local, &T_hi_i](const Alpha& a, const Real& F) {
              // https://en.cppreference.com/w/cpp/language/operator_arithmetic
              const size_t F_i = &F - &Fs.front();
              const auto a_T_hi = Evaluate(F_i, b_s_i_local, T_hi_i);
              return a < a_T_hi;
            }));
    // identify two possible edge cases for CDF
    const auto F_hi_exists = F_hi_i < Fs.size();
    const auto F_lo_exists = F_hi_i != 0;
    if (F_hi_exists && F_lo_exists) {
      // the most common case (2)
      const auto a_F_hi_T_hi = Evaluate(F_hi_i, b_s_i_local, T_hi_i);
      const auto a_F_lo_T_hi = Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i);
      // snap to upper temperature
      const auto a_F_hi = a_F_hi_T_hi;
      const auto a_F_lo = a_F_lo_T_hi;
      // return value of F that would return `a`
      const auto r_a = (a - a_F_lo) / (a_F_hi - a_F_lo);
      const auto F_hi = Fs.at(F_hi_i);
      const auto F_lo = Fs.at(F_hi_i - 1);
      return F_lo + r_a * (F_hi - F_lo);
    }
    else if (F_hi_exists && !F_lo_exists) {
      // use zero alpha at zero CDF (6)
      const auto a_F_hi_T_hi = Evaluate(F_hi_i, b_s_i_local, T_hi_i);
      // snap to upper temperature
      const auto a_F_hi = a_F_hi_T_hi;
      const Alpha a_F_lo = 0;
      // return value of F that would return `a`
      const auto r_a = (a - a_F_lo) / (a_F_hi - a_F_lo);
      const auto F_hi = Fs.at(F_hi_i);
      const Real F_lo = 0;
      return F_lo + r_a * (F_hi - F_lo);
    }
    else if (!F_hi_exists && F_lo_exists) {
      // use alpha cutoff at unity CDF (8)
      const auto a_F_lo_T_hi = Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i);
      // snap to upper temperature
      const auto a_F_hi = alpha_cutoff;
      const auto a_F_lo = a_F_lo_T_hi;
      // return value of F that would return `a`
      const auto r_a = (a - a_F_lo) / (a_F_hi - a_F_lo);
      const Real F_hi = 1;
      const auto F_lo = Fs.at(F_hi_i - 1);
      return F_lo + r_a * (F_hi - F_lo);
    }
    else {
      assert(false);
    }
  }
  else if (!T_hi_exists && T_lo_exists) {
    const size_t F_hi_i = std::distance(
        Fs.cbegin(),
        std::upper_bound(
            Fs.cbegin(), Fs.cend(), a,
            [this, &Fs, &b_s_i_local, &T_hi_i](const Alpha& a, const Real& F) {
              // https://en.cppreference.com/w/cpp/language/operator_arithmetic
              const size_t F_i = &F - &Fs.front();
              const auto a_T_lo = Evaluate(F_i, b_s_i_local, T_hi_i - 1);
              return a < a_T_lo;
            }));
    // identify two possible edge cases for CDF
    const auto F_hi_exists = F_hi_i < Fs.size();
    const auto F_lo_exists = F_hi_i != 0;
    if (F_hi_exists && F_lo_exists) {
      // the most common case (3)
      const auto a_F_hi_T_lo = Evaluate(F_hi_i, b_s_i_local, T_hi_i - 1);
      const auto a_F_lo_T_lo = Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i - 1);
      // snap to lower temperature
      const auto a_F_hi = a_F_hi_T_lo;
      const auto a_F_lo = a_F_lo_T_lo;
      // return value of F that would return `a`
      const auto r_a = (a - a_F_lo) / (a_F_hi - a_F_lo);
      const auto F_hi = Fs.at(F_hi_i);
      const auto F_lo = Fs.at(F_hi_i - 1);
      return F_lo + r_a * (F_hi - F_lo);
    }
    else if (F_hi_exists && !F_lo_exists) {
      // use zero alpha at zero CDF (7)
      const auto a_F_hi_T_lo = Evaluate(F_hi_i, b_s_i_local, T_hi_i - 1);
      // snap to lower temperature
      const auto a_F_hi = a_F_hi_T_lo;
      const Alpha a_F_lo = 0;
      // return value of F that would return `a`
      const auto r_a = (a - a_F_lo) / (a_F_hi - a_F_lo);
      const auto F_hi = Fs.at(F_hi_i);
      const Real F_lo = 0;
      return F_lo + r_a * (F_hi - F_lo);
    }
    else if (!F_hi_exists && F_lo_exists) {
      // use alpha cutoff at unity CDF (9)
      const auto a_F_lo_T_lo = Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i - 1);
      // snap to lower temperature
      const auto a_F_hi = alpha_cutoff;
      const auto a_F_lo = a_F_lo_T_lo;
      // return value of F that would return `a`
      const auto r_a = (a - a_F_lo) / (a_F_hi - a_F_lo);
      const Real F_hi = 1;
      const auto F_lo = Fs.at(F_hi_i - 1);
      return F_lo + r_a * (F_hi - F_lo);
    }
    else {
      assert(false);
    }
  }
  else {
    assert(false);
  }
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
  const Real r = snap_to_min ? 1
                 : snap_to_lower
                     ? 0
                     : (abs_b - betas.at(b_hi_i - 1)) /
                           (betas.at(b_hi_i) - betas.at(b_hi_i - 1));
  const size_t b_s_i =
      r <= std::uniform_real_distribution{}(p.rng) ? b_hi_i - 1 : b_hi_i;

  // compute minimum and maximum possible value of alpha
  const auto sqrt_E = std::sqrt(E);
  const auto awr = target.awr;
  const auto b_sqrt_E_bkT = std::sqrt(E + b * constants::boltzmann * T);
  const auto akT = awr * constants::boltzmann * T;
  // minimum and maximum alpha values at beta
  const Alpha b_a_min = std::pow(sqrt_E - b_sqrt_E_bkT, 2) / akT;
  const Alpha b_a_max = std::pow(sqrt_E + b_sqrt_E_bkT, 2) / akT;

  // Get partition which contains the sampled beta and get local beta index in
  // that partition
  const size_t P_s_i = std::distance(
                           alpha_partition_beta_begins.cbegin(),
                           std::upper_bound(
                               alpha_partition_beta_begins.cbegin(),
                               alpha_partition_beta_begins.cend(), b_s_i)) -
                       1;
  const auto& P_s = alpha_partitions.at(P_s_i);
  const size_t b_s_i_local = b_s_i - alpha_partition_beta_begins.at(P_s_i);
  // CDFs corresponding to b_a_min and b_a_max
  const auto F_min = P_s.FindCDF(b_a_min, b_s_i_local, T, alpha_cutoff);
  const auto F_max = P_s.FindCDF(b_a_max, b_s_i_local, T, alpha_cutoff);
  // sample a CDF value which is scaled to return a result in [b_a_min, b_a_max)
  const auto F =
      F_min + std::uniform_real_distribution{}(p.rng) * (F_max - F_min);
  // find neighboring CDF value indices
  const auto& Fs = P_s.CDF_modes.GetAxis(0);
  const size_t F_hi_i =
      std::distance(Fs.cbegin(), std::upper_bound(Fs.cbegin(), Fs.cend(), F));
  // identify two possible edge cases for CDF
  const auto F_hi_exists = F_hi_i < Fs.size();
  const auto F_lo_exists = F_hi_i != 0;

  // Evaluate nearest temperatures on the sampled beta grid
  const auto& Ts = P_s.beta_T_modes.GetAxis(1);
  // find neighboring temperature value indices
  const size_t T_hi_i =
      std::distance(Ts.cbegin(), std::upper_bound(Ts.cbegin(), Ts.cend(), T));
  // identify two possible edge cases for temperature
  const auto T_hi_exists = T_hi_i < Ts.size();
  const auto T_lo_exists = T_hi_i != 0;

  // there are a total of 9 valid combinations; handle each one
  // TODO: Factor out common terms
  if (F_hi_exists && F_lo_exists && T_hi_exists && T_lo_exists){
    // the most common case (1)
    const auto a_F_hi_T_hi = P_s.Evaluate(F_hi_i, b_s_i_local, T_hi_i);
    const auto a_F_hi_T_lo = P_s.Evaluate(F_hi_i, b_s_i_local, T_hi_i - 1);
    const auto a_F_lo_T_hi = P_s.Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i);
    const auto a_F_lo_T_lo = P_s.Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i - 1);
    // interpolation fractions
    const auto T_hi = Ts.at(T_hi_i);
    const auto T_lo = Ts.at(T_hi_i - 1);
    const auto r_T = (T - T_lo) / (T_hi - T_lo);
    const auto F_hi = Fs.at(F_hi_i);
    const auto F_lo = Fs.at(F_hi_i - 1);
    const auto r_F = (F - F_lo) / (F_hi - F_lo);
    // linearly interpolate in temperature
    const auto a_F_hi = a_F_hi_T_lo + r_T * (a_F_hi_T_hi - a_F_hi_T_lo);
    const auto a_F_lo = a_F_lo_T_lo + r_T * (a_F_lo_T_hi - a_F_lo_T_lo);
    // update indirect effects
    for (auto& indirect_effect : p.indirect_effects) {
      class Visitor : public Perturbation::IndirectEffect::Visitor {
      public:
        // TODO: Optimize variables passed once debugging is done
        Visitor(
            const Nuclide& target, const ThermalScattering::AlphaPartition& P_s,
            const size_t P_s_i, const size_t b_s_i_local, const size_t F_hi_i,
            const size_t T_hi_i, const Real r_T,
            const Real delta_a)
            : target{target}, P_s{P_s}, P_s_i{P_s_i},
              b_s_i_local{b_s_i_local}, F_hi_i{F_hi_i},
              T_hi_i{T_hi_i}, r_T{r_T}, delta_a{delta_a} {}
        void Visit(Perturbation::IndirectEffect::TNSL& indirect_effect)
            const noexcept final {
          // if current TNSL indirect effect is for a Nuclide that does not
          // match the enclosing scope's target Nuclide, skip
          if (&indirect_effect.perturbation.nuclide != &target) {
            return;
          }
          // Compute indirect effects due to each parameter which influences
          // TNSL dynamics. In an attempt to preserve memory locality, indirect
          // effects are computed in the following order:
          //   -# Singular values (diagonal entries of Sigma matrix)
          //   -# CDF coefficients (two rows of U matrix)
          //   -# Temperature coefficients (two rows of V matrix)
          // You should probably profile to see if this looping scheme actually
          // makes a difference.
          const size_t rank = P_s.singular_values.size();
          // compute indirect effects from perturbing singular values
          const size_t partition_offset =
              indirect_effect.perturbation.partition_offsets.at(P_s_i);
          for (size_t r = 0; r < rank; r++) {
            // retrieve coefficients
            const auto u_F_hi = P_s.CDF_modes.at(F_hi_i, r);
            const auto u_F_lo = P_s.CDF_modes.at(F_hi_i - 1, r);
            const auto v_T_hi = P_s.beta_T_modes.at(b_s_i_local, T_hi_i, r);
            const auto v_T_lo = P_s.beta_T_modes.at(b_s_i_local, T_hi_i - 1, r);
            // update indirect effect
            indirect_effect.indirect_effects.at(partition_offset + r) +=
                -ThermalScattering::BilinearForm(
                    -1, +1, u_F_lo * v_T_lo, u_F_lo * v_T_hi, u_F_hi * v_T_lo,
                    u_F_hi * v_T_hi, 1 - r_T, r_T) /
                delta_a;
          }
          // compute indirect effects from perturbing CDF coefficients
          const size_t CDF_offset =
              partition_offset + rank + (F_hi_i - 1) * rank;
          for (size_t r = 0; r < rank; r++){
            // retrieve coefficients
            const auto s_r = P_s.singular_values.at(r);
            const auto v_T_hi = P_s.beta_T_modes.at(b_s_i_local, T_hi_i, r);
            const auto v_T_lo = P_s.beta_T_modes.at(b_s_i_local, T_hi_i - 1, r);
            // update indirect effects
            const auto x = -s_r * (v_T_lo * (1 - r_T) + v_T_hi * r_T) / delta_a;
            indirect_effect.indirect_effects.at(CDF_offset + r) += -x;
            indirect_effect.indirect_effects.at(CDF_offset + rank + r) += x;
          }
          // compute indirect effects from perturbing beta/temperature
          // coefficients
          const auto& Ts = P_s.beta_T_modes.GetAxis(1);
          const size_t beta_T_offset =
              partition_offset + rank + P_s.CDF_modes.size() +
              b_s_i_local * Ts.size() * rank + (T_hi_i - 1) * rank;
          for (size_t r = 0; r < rank; r++){
            // retrieve coefficeints
            const auto s_r = P_s.singular_values.at(r);
            const auto u_F_hi = P_s.CDF_modes.at(F_hi_i, r);
            const auto u_F_lo = P_s.CDF_modes.at(F_hi_i - 1, r);
            // update indirect effects
            const auto x = -s_r * (-u_F_lo + u_F_hi) / delta_a;
            indirect_effect.indirect_effects.at(beta_T_offset + r) +=
                (1 - r_T) * x;
            indirect_effect.indirect_effects.at(beta_T_offset + rank + r) +=
                r_T * x;
          }
        }

      private:
        const Nuclide& target;
        const ThermalScattering::AlphaPartition& P_s;
        const size_t P_s_i;
        const size_t b_s_i_local;
        const size_t F_hi_i;
        const size_t T_hi_i;
        const Real r_T;
        const Real delta_a;
      };
      indirect_effect->Visit(Visitor{
          target, P_s, P_s_i, b_s_i_local, F_hi_i, T_hi_i, r_T,
          a_F_hi - a_F_lo});
    }
    // linearly interpolate in temperature
    return a_F_lo + r_F * (a_F_hi - a_F_lo);
  }
  else if (F_hi_exists && F_lo_exists && T_hi_exists && !T_lo_exists){
    // snap to upper temperature (2)
    const auto a_F_hi_T_hi = P_s.Evaluate(F_hi_i, b_s_i_local, T_hi_i);
    const auto a_F_lo_T_hi = P_s.Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i);
    // interpolation fractions
    const auto F_hi = Fs.at(F_hi_i);
    const auto F_lo = Fs.at(F_hi_i - 1);
    const auto r_F = (F - F_lo) / (F_hi - F_lo);
    // linearly interpolate in CDF
    return a_F_lo_T_hi + r_F * (a_F_hi_T_hi - a_F_lo_T_hi);
  }
  else if (F_hi_exists && F_lo_exists && !T_hi_exists && T_lo_exists){
    // snap to lower temperature (3)
    const auto a_F_hi_T_lo = P_s.Evaluate(F_hi_i, b_s_i_local, T_hi_i - 1);
    const auto a_F_lo_T_lo = P_s.Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i - 1);
    // interpolation fractions
    const auto F_hi = Fs.at(F_hi_i);
    const auto F_lo = Fs.at(F_hi_i - 1);
    const auto r_F = (F - F_lo) / (F_hi - F_lo);
    // linearly interpolate in CDF
    return a_F_lo_T_lo + r_F * (a_F_hi_T_lo - a_F_lo_T_lo);
  }
  else if (F_hi_exists && !F_lo_exists && T_hi_exists && T_lo_exists){
    // use zero alpha at zero CDF (4)
    const auto a_F_hi_T_hi = P_s.Evaluate(F_hi_i, b_s_i_local, T_hi_i);
    const auto a_F_hi_T_lo = P_s.Evaluate(F_hi_i, b_s_i_local, T_hi_i - 1);
    // interpolation fractions
    const auto T_hi = Ts.at(T_hi_i);
    const auto T_lo = Ts.at(T_hi_i - 1);
    const auto r_T = (T - T_lo) / (T_hi - T_lo);
    // linearly interpolate in temperature
    const auto a_F_hi = a_F_hi_T_lo + r_T * (a_F_hi_T_hi - a_F_hi_T_lo);
    // interpolation fraction
    const auto F_hi = Fs.at(F_hi_i);
    const auto r_F = F / F_hi;
    // linearly interpolate in CDF
    return r_F * a_F_hi;
  }
  else if (!F_hi_exists && F_lo_exists && T_hi_exists && T_lo_exists){
    // use alpha cutoff at unity CDF (5)
    const auto a_F_lo_T_hi = P_s.Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i);
    const auto a_F_lo_T_lo = P_s.Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i - 1);
    // interpolation fractions
    const auto T_hi = Ts.at(T_hi_i);
    const auto T_lo = Ts.at(T_hi_i - 1);
    const auto r_T = (T - T_lo) / (T_hi - T_lo);
    // linearly interpolate in temperature
    const auto a_F_hi = alpha_cutoff;
    const auto a_F_lo = a_F_lo_T_lo + r_T * (a_F_lo_T_hi - a_F_lo_T_lo);
    // linearly interpolate in CDF; note F_hi = 1
    const auto F_lo = Fs.at(F_hi_i - 1);
    return a_F_lo + (F - F_lo) / (1 - F_lo) * (a_F_hi - a_F_lo);
  }
  else if (F_hi_exists && !F_lo_exists && T_hi_exists && !T_lo_exists){
    // snap to upper temperature, use zero alpha at zero CDF (6)
    const auto a_F_hi = P_s.Evaluate(F_hi_i, b_s_i_local, T_hi_i);
    // interpolation fraction
    const auto F_hi = Fs.at(F_hi_i);
    const auto r_F = F / F_hi;
    // linearly interpolate in CDF
    return r_F * a_F_hi;
  }
  else if (F_hi_exists && !F_lo_exists && !T_hi_exists && T_lo_exists){
    // snap to lower temperature, use zero alpha at zero CDF (7)
    const auto a_F_hi = P_s.Evaluate(F_hi_i, b_s_i_local, T_hi_i - 1);
    // interpolation fraction
    const auto F_hi = Fs.at(F_hi_i);
    const auto r_F = F / F_hi;
    // linearly interpolate in CDF
    return r_F * a_F_hi;
  }
  else if (!F_hi_exists && F_lo_exists && T_hi_exists && !T_lo_exists){
    // snap to upper temperature, use alpha cutoff at unity CDF (8)
    const auto a_F_hi = alpha_cutoff;
    const auto a_F_lo = P_s.Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i);
    // linearly interpolate in CDF; note F_hi = 1
    const auto F_lo = Fs.at(F_hi_i - 1);
    return a_F_lo + (F - F_lo) / (1 - F_lo) * (a_F_hi - a_F_lo);
  }
  else if (!F_hi_exists && F_lo_exists && !T_hi_exists && T_lo_exists) {
    // snap to lower temperature, use alpha cutoff at unity CDF (9)
    const auto a_F_hi = alpha_cutoff;
    const auto a_F_lo = P_s.Evaluate(F_hi_i - 1, b_s_i_local, T_hi_i - 1);
    // linearly interpolate in CDF; note F_hi = 1
    const auto F_lo = Fs.at(F_hi_i - 1);
    return a_F_lo + (F - F_lo) / (1 - F_lo) * (a_F_hi - a_F_lo);
  }
  else {
    // your nuclear data is bad and you should feel bad
    assert(false);
  }
}
