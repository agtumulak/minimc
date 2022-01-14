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
  /// @details All HDF5 files referenced in the `tsl` node must follow the
  ///          format specified by HDF5DataSet. The MultiIndex levels of each
  ///          HDF5 dataset must correspond to a particular type of value.
  ///
  ///          The total inelastic scattering cross section is compressed using
  ///          proper orthogonal decomposition:
  ///          @f[
  ///            \sigma_\text{inelastic} (E, T)
  ///            = \sum_{n} c_{n} f_{n}(E) g_{n}(T)
  ///          @f]
  ///
  ///          A path to the singular values @f$ c_{n} @f$ are expected under
  ///          the `total_S` attribute of a `tsl` node. The MultiIndex levels
  ///          must be
  ///          1. Expansion order @f$ n @f$
  ///
  ///          A path to the energy dependent functions @f$ f_{n} (E) @f$ are
  ///          expected under the `total_E` attribute. The MultiIndex levels
  ///          must be
  ///          1. Incident neutron energies @f$ E @f$ in MeV. Let the minimum
  ///             and maximum values of @f$ E @f$ provided be @f$
  ///             E_{\text{min}} @f$ and @f$ E_{\text{max}} @f$, respectively.
  ///             The cutoff energy @f$ E_{\text{cutoff}} @f$ will be set to
  ///             @f$ E_\text{max} @f$. If @f$ \sigma_{\text{inelastic}} (E, T)
  ///             @f$ is requested for @f$ E < E_\text{min} @f$, then @f$
  ///             \sigma_\text{inelastic} (E_\text{min}, T ) @f$ will be
  ///             returned.
  ///          2. Expansion order @f$ n @f$
  ///
  ///          A path to the temperature dependent functions @f$ g_{n} (T) @f$
  ///          are expected under the `total_S` attribute. The MultiIndex
  ///          levels must be
  ///          1. Target temperatures @f$ T @f$ in Kelvin. Let the minimum and
  ///             maximum values of @f$ T @f$ provided be @f$ T_\text{min} @f$
  ///             and @f$ T_\text{max} @f$. If @f$ \sigma_{\text{inelastic}}
  ///             (E, T) @f$ is requested for @f$ T < T_\text{min} @f$ or @f$ T
  ///             > T_\text{max} @f$, then @f$ \sigma_{\text{inelastic}} (E,
  ///             T_\text{min}) @f$ or @f$ \sigma_{\text{inelastic}} (E,
  ///             T_\text{max}) @f$ is returned, respectively.
  ///          2. Expansion order @f$ n @f$
  ///
  ///          TODO: update using POD
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
  /// @details Currently implemented using bilinear interpolation in
  ///          energy and temperature.
  /// @todo Avoid calling HDF5DataSet::at for each order since order
  ///       coefficients are adjacent in memory
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept;
  /// @brief The raison d'etre of this class
  void Scatter(Particle& p) const noexcept;

private:
  // Evaluates the inelastic scattering cross section.
  MacroscopicCrossSection
  EvaluateInelastic(const size_t E_index, const size_t T_index) const noexcept;
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
  // Total scattering cross section singular values
  const HDF5DataSet<1> scatter_xs_singular_values;
  // Total scattering cross section energy dependent functions
  const HDF5DataSet<2> scatter_xs_E_dependence;
  // Total scattering cross setion temperature dependent functions
  const HDF5DataSet<2> scatter_xs_T_dependence;
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
  // Atomic weight ratio of target, yes this is a copy rather than giving a
  // pointer to the parent Nuclide
  const Real awr;

public:
  /// @brief Neutrons above the cutoff energy do not get the thermal scattering
  ///        treatment
  const ContinuousEnergy cutoff_energy =
      scatter_xs_E_dependence.GetAxis(0).back();
};
