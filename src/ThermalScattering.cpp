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
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
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
      // IIFE
      alpha_partitions{[&tnsl_node]() {
        std::vector<AlphaPartition> result;
        for (const auto& partition_node : tnsl_node.child("alpha_partitions")) {
          result.emplace_back(partition_node);
        }
        return result;
      }()},
      // IIFE
      betas{[this]() noexcept {
        std::vector<Real> result;
        for (const auto& partition : alpha_partitions) {
          const auto& betas = partition.beta_T_modes.axes.at(0);
          result.insert(result.end(), betas.cbegin(), betas.cend());
        }
        // betas must be monotonically increasing
        Real prev_beta = 0;
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
              result.back() + partition.beta_T_modes.axes.at(0).size());
        }
        assert(result.back() == betas.size());
        return result;
      }()},
      // IIFE
      alpha_partition_offsets{[this]() {
        std::vector<size_t> result{0};
        for (const auto& partition : alpha_partitions) {
          result.emplace_back(
              partition.CDF_modes.size() + partition.singular_values.size() +
              partition.beta_T_modes.size());
        }
        std::partial_sum(result.cbegin(), result.cend(), result.begin());
        return result;
      }()},
      beta_cutoff{tnsl_node.attribute("beta_cutoff").as_double()},
      alpha_cutoff{tnsl_node.attribute("alpha_cutoff").as_double()},
      target{target} {}

size_t ThermalScattering::CountPerturbableParameters() const noexcept {
  return beta_singular_values.size() + beta_CDF_modes.size() +
         beta_E_T_modes.size();
}

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
  // const auto alpha = SampleAlpha(p, beta, E, T);
  // convert to outgoing energy and cosine
  const auto E_p = E + beta * constants::boltzmann * T;
  // const auto mu =
  //     (E + E_p - double{alpha} * target.awr * constants::boltzmann * T) /
  //     (2 * std::sqrt(E * E_p));
  const auto mu = 0.;
  p.Scatter(mu, E_p);
}

//// private

thread_local std::map<const ThermalScattering*, autodiff::VectorXvar>
    ThermalScattering::beta_coeffs = {};

// AlphaPartition

thread_local std::map<
    const ThermalScattering::AlphaPartition*, autodiff::VectorXvar>
    ThermalScattering::AlphaPartition::alpha_coeffs = {};

ThermalScattering::AlphaPartition::AlphaPartition(
    const pugi::xml_node& partition_node)
    : CDF_modes{partition_node.attribute("CDF").as_string()},
      singular_values{partition_node.attribute("S").as_string()},
      beta_T_modes{partition_node.attribute("beta_T").as_string()},
      size{CDF_modes.size() + singular_values.size() + beta_T_modes.size()} {
  // enforce assumptions about CDF values
  const auto& Fs = CDF_modes.axes.at(0);
  if (Fs.front() <= 0) {
    throw std::runtime_error(
        partition_node.path() + ": First CDF (" + std::to_string(Fs.front()) +
        ") must be strictly positive");
  }
  for (size_t i = 1; i < Fs.size(); i++) {
    const auto F_cur = Fs.at(i);
    const auto F_prev = Fs.at(i - 1);
    if (F_prev >= F_cur) {
      throw std::runtime_error(
          partition_node.path() + ": CDF at index " + std::to_string(i) + " (" +
          std::to_string(F_cur) +
          ") must be strictly greater than CDF at index " +
          std::to_string(i - 1) + "(" + std::to_string(F_prev) + ")");
    }
  }
  if (1 <= Fs.back()) {
    throw std::runtime_error(
        partition_node.path() + ": Last CDF (" + std::to_string(Fs.back()) +
        ") must be strictly less than 1");
  }
}

ThermalScattering::Alpha ThermalScattering::AlphaPartition::Evaluate(
    size_t cdf_index, const size_t local_beta_index,
    const size_t T_index) const {
  const size_t CDF_modes_offset =
      singular_values.size() + CDF_modes.GetOffset(cdf_index);
  const size_t beta_T_modes_offset =
      singular_values.size() + CDF_modes.size() +
      beta_T_modes.GetOffset(local_beta_index, T_index);
  Alpha alpha = 0;
  for (size_t order = 0; order < singular_values.axes.at(0).size(); order++) {
    alpha += alpha_coeffs[this][order] *
             alpha_coeffs[this][CDF_modes_offset + order] *
             alpha_coeffs[this][beta_T_modes_offset + order];
  }
  // TODO: Exponentiate
  return alpha;
}

ThermalScattering::Alpha ThermalScattering::AlphaPartition::Sample(
    Particle& p, const Nuclide& target, const size_t partition_offset,
    const size_t b_i, const Real b, const Temperature T,
    const Real a_cutoff) const noexcept {
  // construct thread local members if necessary
  if (alpha_coeffs.find(this) == alpha_coeffs.cend()) {
    alpha_coeffs[this] = autodiff::VectorXvar(size);
    size_t i = 0;
    for (const auto x : singular_values.values) {
      alpha_coeffs[this][i++] = x;
    }
    for (const auto x : CDF_modes.values) {
      alpha_coeffs[this][i++] = x;
    }
    for (const auto x : beta_T_modes.values) {
      alpha_coeffs[this][i++] = x;
    }
  }
  // identify two possible edge cases for temperature
  const auto& Ts = beta_T_modes.axes.at(1);
  const size_t T_upper_i =
      std::distance(Ts.cbegin(), std::upper_bound(Ts.cbegin(), Ts.cend(), T));
  const auto T_hi_i = T_upper_i < Ts.size() ? T_upper_i : T_upper_i - 1;
  const auto T_lo_i = T_upper_i != 0 ? T_upper_i - 1 : T_upper_i;
  const auto T_hi = Ts.at(T_hi_i);
  const auto T_lo = Ts.at(T_lo_i);
  const auto r_T = T_hi_i != T_lo_i ? (T - T_lo) / (T_hi - T_lo) : 0.;
  // create CDF values, appending two endpoints
  const auto Fs = [this]() {
    const auto& original_Fs = CDF_modes.axes.at(0);
    std::vector<CDF> result(original_Fs.size() + 2);
    result.front() = 0;
    for (size_t original_F_i = 0; original_F_i < original_Fs.size();
         original_F_i++) {
      result[original_F_i + 1] = original_Fs[original_F_i];
    }
    result.back() = 1;
    return result;
  }();
  // linearly interpolate alphas in temperature
  const auto [alphas, fs] = [this, &Fs, &b_i, &T_hi_i, &T_lo_i, &r_T,
                             &a_cutoff]() {
    std::vector<Alpha> alphas(Fs.size());
    alphas.front() = 0; // will be set to optimal value later in this function
    for (size_t F_i = 1; F_i < Fs.size() - 1; F_i++) {
      const auto a_T_hi = autodiff::detail::exp(Evaluate(F_i - 1, b_i, T_hi_i));
      const auto a_T_lo = autodiff::detail::exp(Evaluate(F_i - 1, b_i, T_lo_i));
      alphas[F_i] = (1 - r_T) * a_T_lo + r_T * a_T_hi;
    }
    alphas.back() = a_cutoff;
    // compute optimal first value of alpha and corresponding PDFs
    const auto [optimal_a0, fs] = GetOptimalInitial(alphas, Fs);
    alphas[0] = optimal_a0;
    return std::make_tuple(alphas, fs);
  }();
  // evaluate alpha
  const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
  const auto H = p.Sample();
  // compute minimum and maximum allowed values of alpha
  const auto a_min =
      std::pow(std::sqrt(E) - std::sqrt(E + b * constants::boltzmann * T), 2) /
      (target.awr * constants::boltzmann * T);
  const auto a_max =
      std::pow(std::sqrt(E) + std::sqrt(E + b * constants::boltzmann * T), 2) /
      (target.awr * constants::boltzmann * T);
  // compute corresponding CDF values
  const auto F_min = EvaluateQuadratic(alphas, Fs, fs, a_min);
  const auto F_max = EvaluateQuadratic(alphas, Fs, fs, a_max);
  // compute scaled CDF value and return quadratically interpolated alpha
  const auto H_hat = F_min.expr->val + H * (F_max - F_min)->val;
  const auto [alpha, unscaled_f] = SolveQuadratic(alphas, Fs, fs, H_hat);
  const auto f = unscaled_f / (F_max - F_min);
  return alpha;
}

// ThermalScattering

//// private

std::tuple<Real, autodiff::var> ThermalScattering::SolveCubic(
    const std::array<autodiff::var, 4>& xs, const std::array<CDF, 4>& ys,
    const Real y) {
  // 4 neighboring points define 3 intervals. Identify 3 secant lines.
  const auto s0 = (ys[1] - ys[0]) / (xs[1] - xs[0]);
  const auto s1 = (ys[2] - ys[1]) / (xs[2] - xs[1]);
  const auto s2 = (ys[3] - ys[2]) / (xs[3] - xs[2]);
  // Identfy harmonic mean of adjacent secant lines as derivative
  const auto m1 = 2 * s0 * s1 / (s0 + s1);
  const auto m2 = 2 * s1 * s2 / (s1 + s2);
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
  constexpr size_t max_iterations = 20;
  constexpr Real yprevp_eps = 1e-4;
  constexpr Real t_tol = 1e-8;
  Real tprev = 0.5;
  Real tcur;
  for (size_t i = 0; i < max_iterations; i++) {
    const Real yprev = cav * tprev * tprev * tprev + cbv * tprev * tprev +
                   ccv * tprev + (cdv - y);
    const Real yprevp = qav * tprev * tprev + qbv * tprev + qcv;
    if (std::abs(yprevp) < yprevp_eps) {
      throw std::runtime_error("derivative too small");
    }
    tcur = tprev - yprev / yprevp;
    if (std::abs(tcur - tprev) < t_tol) {
      break;
    }
    tprev = tcur;
    if (i == max_iterations - 1) {
      throw std::runtime_error("max iterations exceeded");
    }
  }
  Real beta = tcur * (delta_x->val) + xs[1].expr->val;
  // evaluate derivative of y with respect to x; recall dy(x)/dx = dy(t)/dt
  const auto f = (qa * tcur * tcur + qb * tcur + qc);
  return {beta, f};
}

std::tuple<autodiff::var, std::vector<autodiff::var>>
ThermalScattering::GetOptimalInitial(
    const std::vector<autodiff::var>& xs,
    const std::vector<Real>& ys) noexcept {
  // difference between subsequent values, the first entry is zero
  // because of `std::adjacent_difference` behavior
  const auto delta_xs = [&xs]() {
    std::vector<autodiff::var> result(xs.size());
    std::adjacent_difference(xs.cbegin(), xs.cend(), result.begin());
    return result;
  }();
  // Terms used in the recursion relation used to compute subsequent PDF
  // values. The first value, c0, is set to zero. The second value, c1, is to
  // be determined.
  std::vector<autodiff::var> signed_cs(ys.size());
  signed_cs[0] = 0.; // by definition
  signed_cs[1] = 0.; // to be determined
  for (size_t i = 2; i < signed_cs.size(); i++) {
    signed_cs[i] = std::pow(-1, i) * 2 * (ys[i] - ys[i - 1]) / delta_xs[i];
  }
  // evaluate cummulative sum of `signed_cs`; recall that the first two
  // elements are zero
  const auto cumsum_signed_cs = [&signed_cs]() {
    std::vector<autodiff::var> result(signed_cs.size());
    std::partial_sum(signed_cs.cbegin(), signed_cs.cend(), result.begin());
    return result;
  }();
  // precompute term that appears in both
  const auto delta_xs_pow4 = [&delta_xs]() {
    std::vector<autodiff::var> result(delta_xs.size());
    std::transform(
        delta_xs.cbegin(), delta_xs.cend(), result.begin(),
        [](const auto& delta_xs) {
          return autodiff::detail::pow(delta_xs, 4);
        });
    return result;
  }();
  // compute numerator terms
  const auto numerator_factors = [&cumsum_signed_cs, &signed_cs]() {
    std::vector<autodiff::var> result(cumsum_signed_cs.size());
    std::transform(
        cumsum_signed_cs.cbegin(), cumsum_signed_cs.cend(), signed_cs.cbegin(),
        result.begin(), [](const auto& cumsum_signed_c, const auto& signed_c) {
          return 2 * cumsum_signed_c - signed_c;
        });
    return result;
  }();
  const auto numerator_terms = [&numerator_factors, &delta_xs_pow4]() {
    std::vector<autodiff::var> result(numerator_factors.size());
    std::transform(
        numerator_factors.cbegin(), numerator_factors.cend(),
        delta_xs_pow4.cbegin(), result.begin(),
        std::multiplies<autodiff::var>());
    return result;
  }();
  // assign optimal value of c1
  signed_cs[1] = -0.5 *
                 std::accumulate(
                     std::next(numerator_terms.cbegin(), 2),
                     numerator_terms.cend(), autodiff::var{0.}) /
                 std::accumulate(
                     std::next(delta_xs_pow4.cbegin(), 2), delta_xs_pow4.cend(),
                     autodiff::var{0.});
  // expression for derivatives is similar to `cumsum_signed_cs`
  std::vector<autodiff::var> optimal_derivatives(cumsum_signed_cs);
  for (size_t i = 1; i < cumsum_signed_cs.size(); i++) {
    optimal_derivatives[i] += signed_cs[1];
    optimal_derivatives[i] *= std::pow(-1, i);
  };
  // get optimal value of x0
  const auto optimal_x = xs[1] + 2 * (ys[1] - ys[0]) / signed_cs[1];
  return {optimal_x, optimal_derivatives};
}

autodiff::var ThermalScattering::EvaluateQuadratic(
    const std::vector<autodiff::var>& xs, const std::vector<Real>& ys,
    const std::vector<autodiff::var>& fs, Real x) noexcept {
  // identify first value on x grid which is strictly greater than x
  const size_t hi_i = std::distance(
      xs.cbegin(),
      std::upper_bound(
          xs.cbegin(), xs.cend(), x,
          // this speeds up calculation when autodiff doesn't have to be used
          [](const auto& x, const auto& xs_element) {
            return x < xs_element.expr->val;
          }));
  // if x is strictly less than least value in xs, we say it is zero
  if (hi_i == 0) {
    return 0.;
  }
  // identify quadratic coefficients
  const auto r = (x - xs[hi_i - 1]) / (xs[hi_i] - xs[hi_i - 1]);
  const auto a = 0.5 * (fs[hi_i] - fs[hi_i - 1]) * (xs[hi_i] - xs[hi_i - 1]);
  const auto b = fs[hi_i - 1] * (xs[hi_i] - xs[hi_i - 1]);
  const auto c = ys[hi_i - 1];
  return a * r * r + b * r + c;
}

std::tuple<Real, autodiff::var> ThermalScattering::SolveQuadratic(
    const std::vector<autodiff::var>& xs, const std::vector<Real>& ys,
    const std::vector<autodiff::var>& fs, Real y) noexcept {
  // identify first value on y grid which is strictly greater than y
  const size_t hi_i =
      std::distance(ys.cbegin(), std::upper_bound(ys.cbegin(), ys.cend(), y));
  // identify quadratic coefficients
  const auto delta_x = (xs[hi_i] - xs[hi_i - 1]);
  const auto delta_f = (fs[hi_i] - fs[hi_i - 1]);
  // quadratic coefficients; no need to use autodiff
  const auto a = 0.5 * delta_f->val * delta_x->val;
  const auto b = fs[hi_i - 1].expr->val * delta_x->val;
  const auto c = -(y - ys[hi_i - 1]);
  // solve quadratic equation
  const auto x =
      b >= 0
          ? xs[hi_i - 1].expr->val +
                delta_x->val * (2 * c) / (-b - std::sqrt(b * b - 4 * a * c))
          :
          // handle special case where optimal derivative of first value may be
          // negative
          xs[hi_i - 1].expr->val +
              delta_x->val * (-b + std::sqrt(b * b - 4 * a * c)) / (2 * a);
  // calculate corresponding derivative (may or may not be the PDF just yet)
  // we use non-autodiff value of x since DOS computes the derivative of the
  // probability of sampling the unperturbed value
  const auto f = fs[hi_i - 1] + (x - xs[hi_i - 1]) / delta_x * delta_f;
  return {x, f};
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

ThermalScattering::Beta ThermalScattering::EvaluateBeta(
    const size_t cdf_index, const size_t E_index,
    const size_t T_index) const noexcept {
  const size_t CDF_modes_offset =
      beta_singular_values.size() + beta_CDF_modes.GetOffset(cdf_index);
  const size_t E_T_modes_offset = beta_singular_values.size() +
                                  beta_CDF_modes.size() +
                                  beta_E_T_modes.GetOffset(E_index, T_index);
  Beta beta = 0;
  for (size_t order = 0; order < beta_singular_values.axes.at(0).size();
       order++) {
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
  const auto& Es = beta_E_T_modes.axes.at(0);
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

  // identify two possible edge cases for temperature
  const auto& Ts = beta_E_T_modes.axes.at(1);
  const size_t T_upper_i =
      std::distance(Ts.cbegin(), std::upper_bound(Ts.cbegin(), Ts.cend(), T));
  const auto T_hi_i = T_upper_i < Ts.size() ? T_upper_i : T_upper_i - 1;
  const auto T_lo_i = T_upper_i != 0 ? T_upper_i - 1 : T_upper_i;
  const auto T_hi = Ts.at(T_hi_i);
  const auto T_lo = Ts.at(T_lo_i);
  const auto r_T = T_hi_i != T_lo_i ? (T - T_lo) / (T_hi - T_lo) : 0.;

  const auto& Gs = beta_CDF_modes.axes[0];
  const CDF G = p.Sample();

  // get four neighboring betas in CDF
  const auto [xs, ys] = [&]() {
    std::array<CDF, 4> ys;
    std::array<autodiff::var, 4> xs;
    const size_t hi_i =
        std::distance(Gs.cbegin(), std::upper_bound(Gs.cbegin(), Gs.cend(), G));
    // identifies the CDF index to start reading from dataset
    size_t i_read_begin;
    // identifies one past the last index to read from dataset
    size_t i_read_end;
    // identifies the index to writing interpolated values
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
  const auto [beta, f] = SolveCubic(xs, ys, G);
  // shift result by b_min
  const auto beta_shifted = beta - E / (constants::boltzmann * T);

  // compute sensitivities
  for (auto& indirect_effect : p.indirect_effects) {
    class Visitor : public Perturbation::IndirectEffect::Visitor {
    public:
      Visitor(
          const Nuclide& target, const autodiff::var& f,
          autodiff::VectorXvar& beta_coeffs)
          : target{target}, f{f}, beta_coeffs{beta_coeffs} {}
      void Visit(Perturbation::IndirectEffect::TNSL& indirect_effect)
          const noexcept final {
        // if current TNSL indirect effect is for a Nuclide that does not
        // match the enclosing scope's target Nuclide, skip
        if (&indirect_effect.perturbation.nuclide != &target) {
          return;
        }
        // compute all derivatives in a single reverse pass
        const auto grad = autodiff::gradient(f, beta_coeffs);
        for (auto i = 0; i < grad.size(); i++) {
          indirect_effect.indirect_effects.at(i) += grad[i] / f.expr->val;
        }
      }

    private:
      const Nuclide& target;
      const autodiff::var& f;
      autodiff::VectorXvar& beta_coeffs;
    };

    indirect_effect->Visit(Visitor{target, f, beta_coeffs[this]});
  }
  return beta_shifted;
}

Real ThermalScattering::SampleAlpha(
    Particle& p, const Real& b, ContinuousEnergy E, Temperature T) const {

  // assume S(a,b) = S(a,-b)
  const Real abs_b = std::abs(b);
  // get sign of beta: https://stackoverflow.com/a/4609795/5101335
  const int_fast8_t sgn_b = (0 < b) - (b < 0);
  // identify upper and lower beta gridpoints to use
  const auto [b_hi_i, b_lo_i] = [this, sgn_b, abs_b, E, T]() {
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
  const auto alpha = P_s.Sample(
      p, target, alpha_partition_offsets.at(P_s_i), b_s_i_local, b, T,
      alpha_cutoff);
  return alpha.expr->val;
}
