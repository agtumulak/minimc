#pragma once

#include "BasicTypes.hpp"
#include "ContinuousMap.hpp"
#include "HDF5DataSet.hpp"

#include <cstddef>
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
  // https://en.wikipedia.org/wiki/Bilinear_form#Coordinate_representation
  static Real BilinearForm(
      const Real x0, const Real x1, const Real a00, const Real a01,
      const Real a10, const Real a11, const Real y0, const Real y1) noexcept;
  /// @brief Constructs thermal scattering data from a `tnsl` node and target
  ///        Nuclide
  ThermalScattering(
      const pugi::xml_node& tnsl_node, const Nuclide& target) noexcept;
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
  // Sample an outgoing energy. Requires Particle energy is strictly below
  // ThermalScattering::cutoff_energy. Uses histogram interpolation in PDF.
  Beta SampleBeta(Particle& p, ContinuousEnergy E, Temperature T) const;
  // Sample an outgoing cosine given an outgoing energy.
  Alpha SampleAlpha(
      Particle& p, const Beta& b, ContinuousEnergy E, Temperature T) const;
  // Majorant cross section
  const ContinuousMap<ContinuousEnergy, MicroscopicCrossSection> majorant;
  // Total scattering cross section POD coefficients
  const HDF5DataSet<2> scatter_xs_T;
  const HDF5DataSet<1> scatter_xs_S;
  const HDF5DataSet<2> scatter_xs_E;

public:
  // Encapsulates data required to sample a value of beta using proper
  // orthogonal decomposition; limited to a given incident energy range
  struct BetaPartition {
    // Constructs a BetaPartition from a `partition` node
    BetaPartition(const pugi::xml_node& partition_node);
    // Evaluates beta using proper orthogonal decomposition and linear
    // interpolation in T
    Beta Evaluate(const size_t cdf_index, const size_t E_index, Temperature T)
        const noexcept;
    // Contains CDF modes
    const HDF5DataSet<2> CDF_modes;
    // Contains singular values
    const HDF5DataSet<1> singular_values;
    // Contains energy and temperature modes
    const HDF5DataSet<3> E_T_modes;
  };
  // beta proper orthogonal decomposition coefficients
  const std::vector<BetaPartition> beta_partitions;

private:
  // concatenated vector of all incident energies found in each partition
  const std::vector<ContinuousEnergy> Es;
  // A vector of one-past-the-last energy index for each partition
  const std::vector<size_t> beta_partition_E_ends;

public:
  // Encapsulates data required to sample a value of alpha using proper
  // orthogonal decomposition; limited to a given beta range
  struct AlphaPartition {
    // Constructs an AlphaPartition from a `partition` node
    AlphaPartition(const pugi::xml_node& partition_node);
    // Evaluates alpha using proper orthogonal decomposition and linear
    // interpolation in T
    Alpha Evaluate(
        const size_t cdf_index, const size_t local_beta_index,
        const size_t T_index) const;
    // Returns a CDF value that would return the given value of alpha
    CDF FindCDF(
        const Alpha a, const size_t b_s_i_local, const Temperature T,
        const Alpha alpha_cutoff) const noexcept;
    // Contains CDF modes
    const HDF5DataSet<2> CDF_modes;
    // Contains singular values
    const HDF5DataSet<1> singular_values;
    // Contains beta and temperature modes
    const HDF5DataSet<3> beta_T_modes;
  };
  // alpha proper orthogonal decomposition coefficients
  const std::vector<AlphaPartition> alpha_partitions;

private:
  // concatenated vector of all betas found in each partition_ends
  const std::vector<Beta> betas;
  // Identifies the global beta index where each partition begins
  const std::vector<size_t> alpha_partition_beta_begins;
  // Maximum value of beta which can be sampled
  const Beta beta_cutoff;
  // Maximum value of alpha which can be sampled
  const Alpha alpha_cutoff;
  // the target Nuclide
  const Nuclide& target;

public:
  /// @brief Neutrons above the cutoff energy do not get the thermal scattering
  ///        treatment
  const ContinuousEnergy cutoff_energy = scatter_xs_E.GetAxis(0).back();
};
