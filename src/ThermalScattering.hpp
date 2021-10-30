#pragma once

#include "BasicTypes.hpp"
#include "ContinuousMap.hpp"
#include "HDF5DataSet.hpp"

#include <vector>
#include <cstddef>

namespace pugi {
class xml_node;
}
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
  /// @brief Constructs thermal scattering data from a `tsl` node
  /// @details HDF5 files referenced in the `tsl` node must follow the format
  ///          specified by HDF5DataSet. Additionally, the levels of the
  ///          @f$ \beta @f$ MultiIndex must be in the following order:
  ///          1. Incident neutron energies in @f$(0, E_{\text{cutoff}})@f$
  ///             (non-inclusive) in MeV
  ///          2. CDF values in @f$ (0, 1) @f$ (non-inclusive)
  ///          3. Coefficient labels
  ///
  ///          Also, the levels of the @f$ \alpha @f$ MultiIndex must be in the
  ///          following order:
  ///          1. Sampled energy transfer @f$ \beta @f$ values in @f$ \left[0,
  ///             \max (\frac{E}{kT}, \beta_{\text{cutoff}}) \right) @f$. Note
  ///             that the conditional probability in @f$ \alpha @f$ is
  ///             symmetric about @f$ \beta @f$.
  ///          2. CDF values in @f$ (0, 1) @f$ (non-inclusive).
  ///          3. Coefficient labels
  ///
  ///          The actual values are therefore the coefficients required to
  ///          reconstruct @f$ \beta(T) @f$ or @f$ \alpha(T) @f$.
  ThermalScattering(const pugi::xml_node& tsl_node) noexcept;
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
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept;
  /// @brief The raison d'etre of this class
  void Scatter(Particle& p) const noexcept;

private:
  // Evaluates beta using a functional expansion in temperature
  Real EvaluateBeta(const size_t E_index, const size_t cdf_index, Temperature T)
      const noexcept;
  // Evaluates alpha using a functional expansion in temperature
  Real EvaluateAlpha(
      const size_t beta_index, const size_t cdf_index,
      Temperature T) const noexcept;
  // Sample an outgoing energy. Requires Particle energy is strictly below
  // ThermalScattering::cutoff_energy. Uses histogram interpolation in PDF.
  Beta
  SampleBeta(Particle& p, ContinuousEnergy E, Temperature T) const noexcept;
  // Sample an outgoing cosine given an outgoing energy.
  Alpha SampleAlpha(
      Particle& p, const Beta& b, ContinuousEnergy E,
      Temperature T) const noexcept;
  // Raw cumulative distribution function data for beta and alpha (3
  // independent axes each)
  const HDF5DataSet<3> beta_cdf, alpha_cdf;
  // Majorant cross section
  const ContinuousMap<ContinuousEnergy, MicroscopicCrossSection> majorant;
  // Total scattering cross section
  const ContinuousMap<
      Temperature, ContinuousMap<ContinuousEnergy, MicroscopicCrossSection>>
      total;
  // Incident energies for beta
  const std::vector<ContinuousEnergy> beta_energies = beta_cdf.GetAxis(0);
  // CDF values for beta
  const std::vector<Real> beta_cdfs = beta_cdf.GetAxis(1);
  // Number of beta coefficients
  const size_t n_beta_coefficients = beta_cdf.GetAxis(2).size();
  // Maximum value of beta which can be sampled
  const Beta beta_cutoff;
  // Beta values for energy-independent alpha distributions
  const std::vector<Real> alpha_betas = alpha_cdf.GetAxis(0);
  // CDF valeus for alpha
  const std::vector<Real> alpha_cdfs = alpha_cdf.GetAxis(1);
  // Number of alpha coefficients
  const size_t n_alpha_coefficients = alpha_cdf.GetAxis(2).size();
  // Maximum value of alpha which can be sampled
  const Alpha alpha_cutoff;
  // Range of valid temperatures
  const Temperature min_temperature, max_temperature;
  // Atomic weight ratio of target, yes this is a copy rather than giving a
  // pointer to the parent Nuclide
  const Real awr;

public:
  /// @brief Neutrons above the cutoff energy do not get the thermal scattering
  ///        treatment
  const ContinuousEnergy cutoff_energy = beta_energies.back();
};
