#pragma once

#include "BasicTypes.hpp"
#include "HDF5DataSet.hpp"
#include "pugixml.hpp"

#include <vector>
#include <cstddef>

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
  /// @brief Sample an outgoing energy. Requires Particle energy is strictly
  ///        below ThermalScattering::cutoff_energy. Uses histogram
  ///        interpolation in PDF.
  Beta
  SampleBeta(Particle& p, ContinuousEnergy E, Temperature T) const noexcept;
  /// @brief Sample an outgoing cosine given an outgoing energy.
  Alpha SampleAlpha(
      Particle& p, const Beta& b, ContinuousEnergy E,
      Temperature T) const noexcept;

private:
  // Evaluates beta using a functional expansion in temperature
  Real EvaluateBeta(const size_t E_index, const size_t cdf_index, Temperature T)
      const noexcept;
  // Evaluates alpha using a functional expansion in temperature
  Real EvaluateAlpha(
      const size_t beta_index, const size_t cdf_index,
      Temperature T) const noexcept;
  // Raw cumulative distribution function data for beta and alpha
  const HDF5DataSet beta_cdf, alpha_cdf;
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
  // Atomic weight ratio of target
  const Real awr;

public:
  /// @brief Neutrons above the cutoff energy do not get the thermal scattering
  ///        treatment
  const ContinuousEnergy cutoff_energy = beta_energies.back();
};
