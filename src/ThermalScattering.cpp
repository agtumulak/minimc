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
#include <stdexcept>
#include <type_traits>
#include <variant>

// ThermalScattering

//// public

ThermalScattering::ThermalScattering(
    const pugi::xml_node& tnsl_node, const Nuclide& target)
    : majorant{HDF5DataSet<1>{tnsl_node.attribute("majorant").as_string()}
                   .ToContinuousMap()},
      scatter_xs_T{tnsl_node.attribute("total_T").as_string()},
      scatter_xs_S{tnsl_node.attribute("total_S").as_string()},
      scatter_xs_E{tnsl_node.attribute("total_E").as_string()},
      beta_CDF_modes{tnsl_node.attribute("beta_CDF").as_string()},
      beta_singular_values{tnsl_node.attribute("beta_S").as_string()},
      beta_E_T_modes{tnsl_node.attribute("beta_E_T").as_string()},
      alpha_CDF_modes{tnsl_node.attribute("alpha_CDF").as_string()},
      alpha_singular_values{tnsl_node.attribute("alpha_S").as_string()},
      alpha_beta_T_modes{tnsl_node.attribute("alpha_beta_T").as_string()},
      beta_cutoff{tnsl_node.attribute("beta_cutoff").as_double()},
      alpha_cutoff{tnsl_node.attribute("alpha_cutoff").as_double()},
      target{target} {}

size_t ThermalScattering::CountPerturbableParameters() const noexcept {
  return beta_singular_values.size() + beta_CDF_modes.size() +
         beta_E_T_modes.size() + alpha_singular_values.size() +
         alpha_CDF_modes.size() + alpha_beta_T_modes.size();
}

bool ThermalScattering::IsValid(const Particle& p) const noexcept {
  const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
  return p.GetType() == Particle::Type::neutron && E < cutoff_energy;
}

MicroscopicCrossSection
ThermalScattering::GetCellMajorant(const Particle& p) const noexcept {
  return majorant.at(std::get<ContinuousEnergy>(p.GetEnergy()));
}

MicroscopicCrossSection
ThermalScattering::GetTotal(const Particle& p) const noexcept {

  // Find index of Energy above and below Particle Energy
  const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
  const auto& Es = scatter_xs_E.axes.at(0);
  const size_t E_hi_i =
      std::distance(Es.cbegin(), std::upper_bound(Es.cbegin(), Es.cend(), E));
  // We require that E is strictly less than cutoff energy
  assert(E_hi_i != Es.size());
  // If E < E_min, we use the same index for the lower value of E
  const auto below_E_min = E_hi_i == 0;
  const size_t E_lo_i = below_E_min ? E_hi_i : E_hi_i - 1;

  // Find index of Temperature above and below target Temperature
  const Temperature T = p.GetCell().temperature->at(p.GetPosition());
  const auto& Ts = scatter_xs_T.axes.at(0);
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
  const auto mu =
      (E + E_p - double{alpha} * target.awr * constants::boltzmann * T) /
      (2 * std::sqrt(E * E_p));
  p.Scatter(mu, E_p);
}

//// private

autodiff::var
ThermalScattering::EvaluateG(const autodiff::var& s, const autodiff::var& t) {
  // avoid autodiff whenever possible, generally used when value will not vary
  // continuously
  const auto sv = s.expr->val;
  const auto tv = t.expr->val;
  // adapted from "An Improved Monotone Piecewise Cubic Interpolation
  // Algorithm" by F. N. Fritsch and J. Butland (UCRL 85104)
  constexpr Real m = 3.; // 1 <= m <= 3; 2 is Butland; 3 is recommended
  if (sv * tv > 0) {
    const auto u = std::abs(sv) < std::abs(tv) ? s : t;
    const auto v = std::abs(sv) > std::abs(tv) ? s : t;
    const auto r = u / v;
    const int_fast8_t sgn_s = (0 < sv) - (sv < 0);
    return sgn_s * m * u / (1 + (m - 1) * r);
  }
  else {
    return 0.;
  }
}

autodiff::var ThermalScattering::EvaluateCubic(
    const std::array<autodiff::var, 4>& xs, const std::array<CDF, 4>& ys,
    const Real x) {
  // 4 neighboring points define 3 intervals. Identify 3 secant lines.
  const auto s01 = (ys[1] - ys[0]) / (xs[1] - xs[0]);
  const auto s12 = (ys[2] - ys[1]) / (xs[2] - xs[1]);
  const auto s23 = (ys[3] - ys[2]) / (xs[3] - xs[2]);
  // use neighboring secant lines to compute derivatives at the 2 interior
  // points
  const auto m1 = EvaluateG(s01, s12);
  const auto m2 = EvaluateG(s12, s23);
  // compute cubic coefficients
  const auto delta_x = xs[2] - xs[1];
  const auto delta_y = ys[2] - ys[1];
  const auto ca = -2 * delta_y + delta_x * (m1 + m2);
  const auto cb = 3 * delta_y - delta_x * (2 * m1 + m2);
  const auto cc = delta_x * m1;
  const auto cd = ys[1];
  // evaluate cubic polynomial
  const auto t = (x - xs[1]) / delta_x;
  return ca * t * t * t + cb * t * t + cc * t + cd;
}

std::tuple<Real, autodiff::var> ThermalScattering::SolveCubic(
    const std::array<autodiff::var, 4>& xs, const std::array<CDF, 4>& ys,
    const Real y) noexcept {
  // 4 neighboring points define 3 intervals. Identify 3 secant lines.
  const auto s01 = (ys[1] - ys[0]) / (xs[1] - xs[0]);
  const auto s12 = (ys[2] - ys[1]) / (xs[2] - xs[1]);
  const auto s23 = (ys[3] - ys[2]) / (xs[3] - xs[2]);
  // use neighboring secant lines to compute derivatives at the 2 interior
  // points
  const auto m1 = EvaluateG(s01, s12);
  const auto m2 = EvaluateG(s12, s23);
  // compute cubic coefficients
  const auto delta_x = xs[2] - xs[1];
  const auto delta_y = ys[2] - ys[1];
  const auto ca = -2 * delta_y + delta_x * (m1 + m2);
  const auto cb = 3 * delta_y - delta_x * (2 * m1 + m2);
  const auto cc = delta_x * m1;
  const auto cd = ys[1];
  // compute derivative coefficients
  const auto qa = 3 * ca;
  const auto qb = 2 * cb;
  const auto qc = cc;
  // no need to use autodiff for solving cubic
  const auto cav = (ca)->val; // (c)ubic (a) (v)alue
  const auto cbv = (cb)->val;
  const auto ccv = (cc)->val;
  const auto cdv = cd;
  const auto qav = 3 * cav; // (q)uadratic (a) (v)alue
  const auto qbv = 2 * cbv;
  const auto qcv = ccv;
  // newton raphson method
  // no need to count iterations or check for near-zero derivative because
  // of monotonic interval and bounding of iterations
  constexpr Real nr_epsilon = 1e-12;
  Real tprev = 0.5;
  Real tcur;
  while (true) {
    const Real yprev = cav * tprev * tprev * tprev + cbv * tprev * tprev +
                       ccv * tprev + (cdv - y);
    const Real yprevp = qav * tprev * tprev + qbv * tprev + qcv;
    tcur = tprev - yprev / yprevp;
    // bound next guess to resolve convergence problems for cubics
    tcur = std::max(0., std::min(1., tcur));
    if (std::abs(tcur - tprev) < nr_epsilon) {
      break;
    }
    tprev = tcur;
  }
  Real sampled = tcur * (delta_x->val) + xs[1].expr->val;
  // evaluate derivative of y with respect to x
  const auto t = (sampled - xs[1]) / delta_x;
  const auto p_sampled = (qa * t * t + qb * t + qc) / delta_x;
  return {sampled, p_sampled};
}

MacroscopicCrossSection ThermalScattering::EvaluateInelastic(
    const size_t E_index, const size_t T_index) const noexcept {
  const size_t max_order = scatter_xs_S.axes.at(0).size();
  MicroscopicCrossSection result = 0;
  for (size_t order = 0; order < max_order; order++) {
    result += scatter_xs_S.at(order) * scatter_xs_E.at(E_index, order) *
              scatter_xs_T.at(T_index, order);
  }
  return result;
}

autodiff::var ThermalScattering::EvaluateBeta(
    const size_t cdf_index, const size_t E_index,
    const size_t T_index) const noexcept {
  const size_t CDF_modes_offset =
      beta_singular_values.size() + beta_CDF_modes.GetOffset(cdf_index);
  const size_t E_T_modes_offset = beta_singular_values.size() +
                                  beta_CDF_modes.size() +
                                  beta_E_T_modes.GetOffset(E_index, T_index);
  autodiff::var beta = 0.;
  for (size_t order = 0; order < beta_singular_values.axes[0].size(); order++) {
    beta += beta_coeffs[this][order] *
            beta_coeffs[this][CDF_modes_offset + order] *
            beta_coeffs[this][E_T_modes_offset + order];
  }
  return autodiff::detail::exp(beta);
}

Real ThermalScattering::SampleBeta(
    Particle& p, const ContinuousEnergy E, const Temperature T) const {
  // construct thread local members if necessary
  if (beta_coeffs.find(this) == beta_coeffs.cend()) {
    beta_coeffs[this] = autodiff::VectorXvar(
        beta_singular_values.size() + beta_CDF_modes.size() +
        beta_E_T_modes.size());
    size_t i = 0;
    for (const auto x : beta_singular_values.values) {
      beta_coeffs[this][i++] = x;
    }
    for (const auto x : beta_CDF_modes.values) {
      beta_coeffs[this][i++] = x;
    }
    for (const auto x : beta_E_T_modes.values) {
      beta_coeffs[this][i++] = x;
    }
  }
  // find index of energy value strictly greater than E
  const auto& Es = beta_E_T_modes.axes[0];
  const size_t E_hi_i =
      std::distance(Es.cbegin(), std::upper_bound(Es.cbegin(), Es.cend(), E));
  // We require that E is strictly less than cutoff energy
  if (E_hi_i == Es.size()) {
    throw std::runtime_error(
        std::string{"Requested energy ("} + std::to_string(E) +
        ") must be strictly less than cutoff energy (" +
        std::to_string(cutoff_energy) + ")");
  }
  // Determine interpolation factor. If E < E_min, we force the use of the grid
  // at E_min by setting r = 1. This way, index of the sampled energy E_s_i is
  // always valid.
  const Real r = E_hi_i != 0 ? (E - Es.at(E_hi_i - 1)) /
                                   (Es.at(E_hi_i) - Es.at(E_hi_i - 1))
                             : 1;
  const size_t E_s_i =
      r <= std::uniform_real_distribution{}(p.rng) ? E_hi_i - 1 : E_hi_i;
  // identify two possible edge cases for temperature
  const auto& Ts = beta_E_T_modes.axes[1];
  const size_t T_upper_i =
      std::distance(Ts.cbegin(), std::upper_bound(Ts.cbegin(), Ts.cend(), T));
  const auto T_hi_i = T_upper_i < Ts.size() ? T_upper_i : T_upper_i - 1;
  const auto T_lo_i = T_upper_i != 0 ? T_upper_i - 1 : T_upper_i;
  const auto T_hi = Ts.at(T_hi_i);
  const auto T_lo = Ts.at(T_lo_i);
  const auto r_T = T_hi_i != T_lo_i ? (T - T_lo) / (T_hi - T_lo) : 0.;
  // sample a CDF value
  const auto& Gs = beta_CDF_modes.axes[0];
  const CDF G = p.Sample();
  // get four neighboring (temperature interpolated) betas and corresponding
  // CDFs
  const auto [xs, ys] = [&]() {
    std::array<CDF, 4> ys;
    std::array<autodiff::var, 4> xs;
    const size_t hi_i =
        std::distance(Gs.cbegin(), std::upper_bound(Gs.cbegin(), Gs.cend(), G));
    // CDF/beta index to start reading from dataset
    size_t i_read_begin;
    // one past the last index to read from dataset
    size_t i_read_end;
    // index to start writing tempearture interpolated values
    size_t i_write_begin;
    // handle edge cases near beginning of dataset
    if (hi_i == 0) {
      // mirror reflect first nonzero point
      const auto first_x = (1 - r_T) * EvaluateBeta(0, E_s_i, T_lo_i) +
                           r_T * EvaluateBeta(0, E_s_i, T_hi_i);
      xs[0] = -first_x;
      ys[0] = -Gs.front();
      // beta data is offset from b_min; beta at zero CDF must be zero
      xs[1] = 0.;
      ys[1] = 0.;
      // assign indices
      i_read_begin = 0;
      i_write_begin = 2;
    }
    else if (hi_i == 1) {
      // beta data is offset from b_min; beta at zero CDF must be zero
      xs[0] = 0.;
      ys[0] = 0.;
      // assign indices
      i_read_begin = 0;
      i_write_begin = 1;
    }
    else {
      i_read_begin = hi_i - 2;
      i_write_begin = 0;
    }
    // handle edge cases near end of dataset
    if (hi_i == Gs.size() - 1) {
      // beta data is offset from b_min; beta at unity CDF must equal cutoff
      // after adding b_min
      xs[3] = beta_cutoff + E / (constants::boltzmann * T);
      ys[3] = 1.;
      // assign indices
      i_read_end = Gs.size();
    }
    else if (hi_i == Gs.size()) {
      // beta data is offset from b_min; beta at unity CDF must equal cutoff
      // after adding b_min
      xs[2] = beta_cutoff + E / (constants::boltzmann * T);
      ys[2] = 1.;
      // mirror reflect last point
      const auto last_x =
          (1 - r_T) * EvaluateBeta(Gs.size() - 1, E_s_i, T_lo_i) +
          r_T * EvaluateBeta(Gs.size() - 1, E_s_i, T_hi_i);
      xs[3] = xs[2] + (xs[2] - last_x);
      ys[3] = 1 + (1 - Gs.back());
      // assign indices
      i_read_end = Gs.size();
    }
    else {
      i_read_end = hi_i + 2;
    }
    // write elements to array
    for (auto [i_read, i_write] = std::tuple{i_read_begin, i_write_begin};
         i_read < i_read_end; i_read++, i_write++) {
      xs[i_write] = (1 - r_T) * EvaluateBeta(i_read, E_s_i, T_lo_i) +
                    r_T * EvaluateBeta(i_read, E_s_i, T_hi_i);
      ys[i_write] = Gs[i_read];
    }
    return std::make_tuple(xs, ys);
  }();
  // perform monotone cubic interpolation
  const auto [beta, p_beta] = SolveCubic(xs, ys, G);
  // shift result by b_min
  const auto beta_shifted = beta - E / (constants::boltzmann * T);
  // compute sensitivities
  for (auto& indirect_effect : p.indirect_effects) {
    class Visitor : public Perturbation::IndirectEffect::Visitor {
    public:
      Visitor(
          const Nuclide& target, const autodiff::var& p_beta,
          autodiff::VectorXvar& beta_coeffs)
          : target{target}, p_beta{p_beta}, beta_coeffs{beta_coeffs} {}
      void Visit(Perturbation::IndirectEffect::TNSL& indirect_effect)
          const noexcept final {
        // if current TNSL indirect effect is for a Nuclide that does not
        // match the enclosing scope's target Nuclide, skip
        if (&indirect_effect.perturbation.nuclide != &target) {
          return;
        }
        // compute all derivatives in a single reverse pass
        const auto grad = autodiff::gradient(p_beta, beta_coeffs);
        for (auto i = 0; i < grad.size(); i++) {
          indirect_effect.indirect_effects.at(i) += grad[i] / p_beta.expr->val;
        }
      }

    private:
      const Nuclide& target;
      const autodiff::var& p_beta;
      autodiff::VectorXvar& beta_coeffs;
    };

    indirect_effect->Visit(Visitor{target, p_beta, beta_coeffs[this]});
  }
  return beta_shifted;
}

template <typename T>
T ThermalScattering::EvaluateAlpha(
    size_t cdf_index, const size_t beta_index,
    const size_t T_index) const noexcept {
  static_assert(
      std::is_same_v<T, autodiff::var> || std::is_same_v<T, Real>,
      "Unsupported template type");
  if constexpr (std::is_same_v<T, autodiff::var>) {
    const size_t CDF_modes_offset =
        alpha_singular_values.size() + alpha_CDF_modes.GetOffset(cdf_index);
    const size_t beta_T_modes_offset =
        alpha_singular_values.size() + alpha_CDF_modes.size() +
        alpha_beta_T_modes.GetOffset(beta_index, T_index);
    autodiff::var alpha = 0.;
    for (size_t order = 0; order < alpha_singular_values.axes[0].size();
         order++) {
      alpha += alpha_coeffs[this][order] *
               alpha_coeffs[this][CDF_modes_offset + order] *
               alpha_coeffs[this][beta_T_modes_offset + order];
    }
    return autodiff::detail::exp(alpha);
  }
  else if constexpr (std::is_same_v<T, Real>) {
    Alpha alpha = 0.;
    const size_t CDF_modes_offset = alpha_CDF_modes.GetOffset(cdf_index);
    const size_t beta_T_modes_offset =
        alpha_beta_T_modes.GetOffset(beta_index, T_index);
    for (size_t order = 0; order < alpha_singular_values.axes[0].size();
         order++) {
      alpha += alpha_singular_values.values[order] *
               alpha_CDF_modes.values[CDF_modes_offset + order] *
               alpha_beta_T_modes.values[beta_T_modes_offset + order];
    }
    return std::exp(alpha);
  }
}

autodiff::var ThermalScattering::FindAlphaCDF(
    const Alpha a, const size_t b_s_i, const size_t T_hi_i, const size_t T_lo_i,
    const Real r_T) const noexcept {
  // find first value on alpha grid strictly greater than a
  const auto& Hhats = alpha_CDF_modes.axes[0];
  const auto hi_i = std::distance(
      Hhats.cbegin(),
      std::upper_bound(
          Hhats.cbegin(), Hhats.cend(), a,
          [this, &Hhats, &b_s_i, &T_hi_i, &T_lo_i,
           &r_T](const auto& a, const auto& Hhat) {
            // https://en.cppreference.com/w/cpp/language/operator_arithmetic
            const size_t cdf_i = &Hhat - &Hhats.front();
            // perturbations will not affect `hi_i`, autodiff not needed
            const auto a_T_hi = EvaluateAlpha<Real>(cdf_i, b_s_i, T_hi_i);
            const auto a_T_lo = EvaluateAlpha<Real>(cdf_i, b_s_i, T_lo_i);
            const auto a_T = (1 - r_T) * a_T_lo + r_T * a_T_hi;
            return a < a_T;
          }));
  const auto [xs, ys] = GetNeighboringAlphas(hi_i, b_s_i, T_hi_i, T_lo_i, r_T);
  return EvaluateCubic(xs, ys, a);
}

std::tuple<std::array<autodiff::var, 4>, std::array<CDF, 4>>
ThermalScattering::GetNeighboringAlphas(
    const size_t hi_i, const size_t b_s_i, const size_t T_hi_i,
    const size_t T_lo_i, const Real r_T) const noexcept {
  const auto& Hhats = alpha_CDF_modes.axes[0];
  std::array<CDF, 4> ys;
  std::array<autodiff::var, 4> xs;
  // CDF/alpha index to start reading from dataset
  size_t i_read_begin;
  // one past the last index to read from dataset
  size_t i_read_end;
  // index to start writing tempearture interpolated values
  size_t i_write_begin;
  // handle edge cases near beginning of dataset
  if (hi_i == 0) {
    // mirror reflect first nonzero point
    const auto first_x =
        (1 - r_T) * EvaluateAlpha<autodiff::var>(0, b_s_i, T_lo_i) +
        r_T * EvaluateAlpha<autodiff::var>(0, b_s_i, T_hi_i);
    xs[0] = -first_x;
    ys[0] = -Hhats.front();
    // alpha at zero CDF is defined to be zero
    xs[1] = 0.;
    ys[1] = 0.;
    // assign indices
    i_read_begin = 0;
    i_write_begin = 2;
  }
  else if (hi_i == 1) {
    // alpha at zero CDF is defined to be zero
    xs[0] = 0.;
    ys[0] = 0.;
    // assign indices
    i_read_begin = 0;
    i_write_begin = 1;
  }
  else {
    i_read_begin = hi_i - 2;
    i_write_begin = 0;
  }
  // handle edge cases near end of dataset
  if (hi_i == Hhats.size() - 1) {
    // alpha at unity CDF is defiend to be alpha_cutoff
    xs[3] = alpha_cutoff;
    ys[3] = 1.;
    // assign indices
    i_read_end = Hhats.size();
  }
  else if (hi_i == Hhats.size()) {
    // alpha at unity CDF is defiend to be alpha_cutoff
    xs[2] = alpha_cutoff;
    ys[2] = 1.;
    // mirror reflect last point
    const auto last_x =
        (1 - r_T) *
            EvaluateAlpha<autodiff::var>(Hhats.size() - 1, b_s_i, T_lo_i) +
        r_T * EvaluateAlpha<autodiff::var>(Hhats.size() - 1, b_s_i, T_hi_i);
    xs[3] = xs[2] + (xs[2] - last_x);
    ys[3] = 1 + (1 - Hhats.back());
    // assign indices
    i_read_end = Hhats.size();
  }
  else {
    i_read_end = hi_i + 2;
  }
  // write elements to array
  for (auto [i_read, i_write] = std::tuple{i_read_begin, i_write_begin};
       i_read < i_read_end; i_read++, i_write++) {
    xs[i_write] =
        (1 - r_T) * EvaluateAlpha<autodiff::var>(i_read, b_s_i, T_lo_i) +
        r_T * EvaluateAlpha<autodiff::var>(i_read, b_s_i, T_hi_i);
    ys[i_write] = Hhats[i_read];
  }
  return {xs, ys};
}

Real ThermalScattering::SampleAlpha(
    Particle& p, const Real& b, ContinuousEnergy E, Temperature T) const {
  // construct thread local members if necessary
  if (alpha_coeffs.find(this) == alpha_coeffs.cend()) {
    alpha_coeffs[this] = autodiff::VectorXvar(
        alpha_singular_values.size() + alpha_CDF_modes.size() +
        alpha_beta_T_modes.size());
    size_t i = 0;
    for (const auto x : alpha_singular_values.values) {
      alpha_coeffs[this][i++] = x;
    }
    for (const auto x : alpha_CDF_modes.values) {
      alpha_coeffs[this][i++] = x;
    }
    for (const auto x : alpha_beta_T_modes.values) {
      alpha_coeffs[this][i++] = x;
    }
  }
  // assume S(a,b) = S(a,-b)
  const Real abs_b = std::abs(b);
  // get sign of beta: https://stackoverflow.com/a/4609795/5101335
  const int_fast8_t sgn_b = (0 < b) - (b < 0);
  // identify upper and lower beta gridpoints to use
  const auto& betas = alpha_beta_T_modes.axes.at(0);
  const auto [b_hi_i, b_lo_i] = [this, sgn_b, abs_b, betas, E, T]() {
    // find index of beta value strictly greater than abs_b
    const size_t b_upper_i = std::distance(
        betas.cbegin(), std::upper_bound(betas.cbegin(), betas.cend(), abs_b));
    // if only one datapoint, we are forced to use it for both cases
    if (betas.size() == 1) {
      return std::make_tuple(size_t{0}, size_t{0});
    }
    // requested beta is greater than any available beta; snap to last beta
    if (b_upper_i == betas.size()) {
      return std::make_tuple(b_upper_i - 1, b_upper_i - 1);
    }
    // requested beta is less than any available beta; snap to first value
    if (b_upper_i == 0) {
      return std::make_tuple(size_t{0}, size_t{0});
    }
    // upper beta exists, but it corresponds to physically impossible value
    // so snap to lower value
    const auto physical_upper =
        sgn_b == 1 ? beta_cutoff : E / (constants::boltzmann * T);
    if (physical_upper < betas.at(b_upper_i)) {
      return std::make_tuple(b_upper_i - 1, b_upper_i - 1);
    }
    return std::make_tuple(b_upper_i, b_upper_i - 1);
  }();
  // Determine interpolation factor and sample a beta grid to use, b_s
  const Real r = b_hi_i == b_lo_i ? 0
                                  : (abs_b - betas.at(b_lo_i)) /
                                        (betas.at(b_hi_i) - betas.at(b_lo_i));
  const size_t b_s_i =
      r <= std::uniform_real_distribution{}(p.rng) ? b_lo_i : b_hi_i;
  // identify two possible edge cases for temperature
  const auto& Ts = alpha_beta_T_modes.axes[1];
  const size_t T_upper_i =
      std::distance(Ts.cbegin(), std::upper_bound(Ts.cbegin(), Ts.cend(), T));
  const auto T_hi_i = T_upper_i < Ts.size() ? T_upper_i : T_upper_i - 1;
  const auto T_lo_i = T_upper_i != 0 ? T_upper_i - 1 : T_upper_i;
  const auto T_hi = Ts.at(T_hi_i);
  const auto T_lo = Ts.at(T_lo_i);
  const auto r_T = T_hi_i != T_lo_i ? (T - T_lo) / (T_hi - T_lo) : 0.;
  // identify CDFs corresponding to a_min and a_max
  const auto a_min =
      std::pow(std::sqrt(E) - std::sqrt(E + b * constants::boltzmann * T), 2) /
      (target.awr * constants::boltzmann * T);
  const auto a_max =
      std::pow(std::sqrt(E) + std::sqrt(E + b * constants::boltzmann * T), 2) /
      (target.awr * constants::boltzmann * T);
  const auto Hhat_min = FindAlphaCDF(a_min, b_s_i, T_hi_i, T_lo_i, r_T);
  const auto Hhat_max = FindAlphaCDF(a_max, b_s_i, T_hi_i, T_lo_i, r_T);
  // sample a scaled CDF value
  const auto H = p.Sample();
  const auto Hhat =
      Hhat_min.expr->val + H * (Hhat_max.expr->val - Hhat_min.expr->val);
  const auto& Hhats = alpha_CDF_modes.axes[0];
  const size_t Hhat_hi_i = std::distance(
      Hhats.cbegin(), std::upper_bound(Hhats.cbegin(), Hhats.cend(), Hhat));
  const auto [xs, ys] =
      GetNeighboringAlphas(Hhat_hi_i, b_s_i, T_hi_i, T_lo_i, r_T);
  // perform monotone cubic interpolation
  const auto [alpha, unscaled_p_alpha] = SolveCubic(xs, ys, Hhat);
  const auto p_alpha = unscaled_p_alpha / (Hhat_max - Hhat_min);
  // compute sensitivities
  for (auto& indirect_effect : p.indirect_effects) {
    class Visitor : public Perturbation::IndirectEffect::Visitor {
    public:
      Visitor(
          const Nuclide& target, const autodiff::var& p_alpha,
          autodiff::VectorXvar& alpha_coeffs, const size_t offset)
          : target{target}, p_alpha{p_alpha}, alpha_coeffs{alpha_coeffs},
            offset{offset} {}
      void Visit(Perturbation::IndirectEffect::TNSL& indirect_effect)
          const noexcept final {
        // if current TNSL indirect effect is for a Nuclide that does not
        // match the enclosing scope's target Nuclide, skip
        if (&indirect_effect.perturbation.nuclide != &target) {
          return;
        }
        // compute all derivatives in a single reverse pass
        const auto grad = autodiff::gradient(p_alpha, alpha_coeffs);
        for (auto i = 0; i < grad.size(); i++) {
          indirect_effect.indirect_effects.at(offset + i) +=
              grad[i] / p_alpha.expr->val;
        }
      }

    private:
      const Nuclide& target;
      const autodiff::var& p_alpha;
      autodiff::VectorXvar& alpha_coeffs;
      const size_t offset;
    };

    indirect_effect->Visit(Visitor{
        target, p_alpha, alpha_coeffs[this],
        beta_singular_values.size() + beta_CDF_modes.size() +
            beta_E_T_modes.size()});
  }
  return alpha;
}

// initialize static members

thread_local std::map<const ThermalScattering*, autodiff::VectorXvar>
    ThermalScattering::beta_coeffs = {};

thread_local std::map<const ThermalScattering*, autodiff::VectorXvar>
    ThermalScattering::alpha_coeffs = {};
