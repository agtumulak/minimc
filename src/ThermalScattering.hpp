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
  ///          proper orthogonal decomposition
  ///          @f[
  ///            \sigma_\text{inelastic} (E, T)
  ///            = \sum_{n} c_{n} f_{n}(E) g_{n}(T)
  ///          @f]
  ///          where @f$ E @f$ is the incident neutron energy, and @f$ T @f$ is
  ///          the target temperature.
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
  ///          The sampled value of @f$ \beta @f$ is compressed using proper
  ///          orthogonal decomposition:
  ///          @f[
  ///            \beta (E, F, T)
  ///            = \sum_{n} c_{n} f_{n}(E, F) g_{n}(T)
  ///          @f]
  ///          where @f$ E @f$ is the incident neutron energy, @f$ F @f$ is the
  ///          CDF value, and @f$ T @f$ is the target temperature.
  ///
  ///          A path to the singular values @f$ c_{n} @f$ are expected under
  ///          the `beta_S` attribute of a `tsl` node. The MultiIndex levels
  ///          must be
  ///          1. Expansion order @f$ n @f$
  ///
  ///          A path to the energy and CDF dependent functions @f$ f_{n}
  ///          (E, F) @f$ are expected under the `beta_E_CDF` attribute. The
  ///          MultiIndex levels must be
  ///          1. Incident neutron energies @f$ E @f$ in MeV. Does not need to
  ///             be identical to other energy grids.
  ///          2. CDF values in @f$ (0, 1) @f$ (non-inclusive)
  ///          3. Expansion order @f$ n @f$
  ///
  ///          A path to the temperature dependent functions @f$ g_{n} (T) @f$
  ///          are expected under the `beta_T` attribute. The MultiIndex levels
  ///          must be
  ///          1. Target temperatures @f$ T @f$ in Kelvin. Does not need to be
  ///             identical to other temperature grids.
  ///          2. Expansion order @f$ n @f$
  ///
  ///          The sampled value @f$ \alpha @f$ is compressed using proper
  ///          orthogonal decomposition:
  ///          @f[
  ///            \alpha (\beta, F, T)
  ///            = \sum_{n} c_{n} f_{n}(\beta, F) g_{n}(T)
  ///          @f]
  ///          where @f$ \beta @f$ is the sampled value of @f$ \beta @f$, @f$ F
  ///          @f$ is the CDF value, and @f$ T @f$ is the target temperature.
  ///
  ///          A path to the singular values @f$ c_{n} @f$ are expected under
  ///          the `alpha_S` attribute of a `tsl` node. The MultiIndex levels
  ///          must be
  ///          1. Expansion order @f$ n @f$
  ///
  ///          A path to the beta and CDF dependent functions @f$ f_{n}
  ///          (\beta, F) @f$ are expected under the `alpha_beta_CDF` attribute.
  ///          The MultiIndex levels must be
  ///          1. Sampled beta values @f$ \beta @f$. Does not need to be
  ///             identical to other beta grids.
  ///          2. CDF values in @f$ (0, 1) @f$ (non-inclusive)
  ///          3. Expansion order @f$ n @f$
  ///
  ///          A path to the temperature dependent functions @f$ g_{n} (T) @f$
  ///          are expected under the `alpha_T` attribute. The MultiIndex
  ///          levels must be
  ///          1. Target temperatures @f$ T @f$ in Kelvin. Does not need to be
  ///             identical to other temperature grids.
  ///          2. Expansion order @f$ n @f$
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
  // Evaluates beta using proper orthogonal decomposition and linear
  // interpolation in T
  Beta EvaluateBeta(const size_t E_index, const size_t cdf_index, Temperature T)
      const noexcept;
  // Evaluates alpha using proper orthogonal decomposition and linear
  // interpolation in T
  Alpha EvaluateAlpha(
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
  // Majorant cross section
  const ContinuousMap<ContinuousEnergy, MicroscopicCrossSection> majorant;
  // Total scattering cross section POD coefficients
  const HDF5DataSet<2> scatter_xs_T;
  const HDF5DataSet<1> scatter_xs_S;
  const HDF5DataSet<2> scatter_xs_E;
  // beta proper orthogonal decomposition coefficients
  const HDF5DataSet<2> beta_T;
  const HDF5DataSet<1> beta_S;
  const HDF5DataSet<3> beta_E_CDF;
  // alpha proper orthogonal decomposition coefficients
  const HDF5DataSet<2> alpha_T;
  const HDF5DataSet<1> alpha_S;
  const HDF5DataSet<3> alpha_beta_CDF;
  // Maximum value of beta which can be sampled
  const Beta beta_cutoff;
  // Maximum value of alpha which can be sampled
  const Alpha alpha_cutoff;
  // Atomic weight ratio of target, yes this is a copy rather than giving a
  // pointer to the parent Nuclide
  const Real awr;

public:
  /// @brief Neutrons above the cutoff energy do not get the thermal scattering
  ///        treatment
  const ContinuousEnergy cutoff_energy =
      scatter_xs_E.GetAxis(0).back();
};
