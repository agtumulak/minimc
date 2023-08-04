#pragma once

#include "BasicTypes.hpp"
#include "ContinuousMap.hpp"
#include "HDF5DataSet.hpp"
#include "autodiff/reverse/var/eigen.hpp"
#include "autodiff/reverse/var/var.hpp"

#include <array>
#include <cstddef>
#include <map>
#include <tuple>
#include <vector>

namespace pugi {
class xml_node;
}
class Nuclide;
class Particle;

/// @brief Models temperature-dependent thermal scattering data S(a,b,T)
/// @details Based on
///          <a href="https://doi.org/10.1016/j.anucene.2014.04.028">
///          On-the-fly sampling of temperature-dependent thermal neutron
///          scattering data for Monte Carlo simulations</a>
///          by Pavlou and Ji. For a given energy and CDF value, the value of
///          beta is given as a functional expansion in temperature. For a
///          given beta and CDF value, the value of alpha is given as a
///          functional expansion in temperature.
class ThermalScattering {
public:
  /// @brief Dimensionless energy transfer
  using Beta = Real;
  /// @brief Dimensionless momentum transfer
  using Alpha = Real;
  /// @brief Constructs thermal scattering data from a `tnsl` node and target
  ///        Nuclide
  ThermalScattering(const pugi::xml_node& tnsl_node, const Nuclide& target);
  /// @brief Returns the number of perturbable parameters
  size_t CountPerturbableParameters() const noexcept;
  /// @brief Returns true if Particle is Type::neutron and is strictly below
  ///        the cutoff energy
  bool IsValid(const Particle& p) const noexcept;
  /// @brief Returns the majorant cross section
  /// @details For thermal inelastic scattering, this is larger than @f$
  ///          \Sigma_{\text{elastic}}(E) + \Sigma_{\text{inelastic}}(E) @f$
  ///          because @f$ \beta = 0 @f$ corresponds to elastic scattering.
  /// @todo Throw exception when temperature is out of range
  MicroscopicCrossSection GetMajorant(const Particle& p) const noexcept;
  /// @brief Returns the total cross section
  /// @details Currently implemented using bilinear interpolation in
  ///          energy and temperature.
  /// @todo Avoid calling HDF5DataSet::at for each order since order
  ///       coefficients are adjacent in memory
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept;
  /// @brief The raison d'etre of this class
  void Scatter(Particle& p) const noexcept;

private:
  // compute derivative at point using adjacent secant lines using method by
  // F. N. Fritsch and J. Butland
  static autodiff::var
  EvaluateG(const autodiff::var& s, const autodiff::var& t);
  // evaluate cubic Hermite spline using J. Butland's method
  static autodiff::var EvaluateCubic(
      const std::array<autodiff::var, 4>& xs, const std::array<CDF, 4>& ys,
      const Real x);
  // sample cubic Hermite spline using J. Butland's method
  static std::tuple<Real, autodiff::var> SolveCubic(
      const std::array<autodiff::var, 4>& xs, const std::array<CDF, 4>& ys,
      const Real y) noexcept;
  // evaluates the inelastic scattering cross section.
  MacroscopicCrossSection
  EvaluateInelastic(const size_t E_index, const size_t T_index) const noexcept;
  // evaluates beta using proper orthogonal decomposition
  autodiff::var EvaluateBeta(
      const size_t cdf_index, const size_t E_index,
      const size_t T_index) const noexcept;
  // samples an outgoing energy
  Real SampleBeta(Particle& p, ContinuousEnergy E, Temperature T) const;
  // evaluates alpha using proper orthogonal decomposition
  template <typename T>
  T EvaluateAlpha(
      const size_t cdf_index, const size_t beta_index,
      const size_t T_index) const noexcept;
  // search dataset for CDF value that returns requested alpha
  autodiff::var FindAlphaCDF(
      const Alpha a, const size_t b_s_i, const size_t T_hi_i,
      const size_t T_lo_i, const Real r_T) const noexcept;
  // returns four neighboring alpha values suitable for cubic interpolation
  std::tuple<std::array<autodiff::var, 4>, std::array<CDF, 4>>
  GetNeighboringAlphas(
      const size_t hi_i, const size_t b_s_i, const size_t T_hi_i,
      const size_t T_lo_i, const Real r_T) const noexcept;
  // samples an outgoing cosine given an outgoing energy
  Real SampleAlpha(
      Particle& p, const Real& b, ContinuousEnergy E, Temperature T) const;
  // majorant cross section
  const ContinuousMap<ContinuousEnergy, MicroscopicCrossSection> majorant;
  // total scattering cross section POD coefficients
  const HDF5DataSet<2> scatter_xs_T;
  const HDF5DataSet<1> scatter_xs_S;
  const HDF5DataSet<2> scatter_xs_E;
  // beta POD coefficients; expanded elements are log of offset from minimum
  // value at given energy and temperature
  const HDF5DataSet<2> beta_CDF_modes;
  const HDF5DataSet<1> beta_singular_values;
  const HDF5DataSet<3> beta_E_T_modes;
  // alpha POD coefficients, expanded elements are log of value
  const HDF5DataSet<2> alpha_CDF_modes;
  const HDF5DataSet<1> alpha_singular_values;
  const HDF5DataSet<3> alpha_beta_T_modes;
  // maximum value of beta which can be sampled
  const Real beta_cutoff;
  // maximum value of alpha which can be sampled
  const Real alpha_cutoff;
  // the target Nuclide
  const Nuclide& target;
  // neutrons above the cutoff energy do not get thermal scattering treatment
  const ContinuousEnergy cutoff_energy = scatter_xs_E.axes.at(0).back();
  // each thread needs its own leaf nodes in expression trees
  thread_local static std::map<const ThermalScattering*, autodiff::VectorXvar>
      beta_coeffs, alpha_coeffs;
};
