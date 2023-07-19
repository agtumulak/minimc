#pragma once

#include "BasicTypes.hpp"
#include "ContinuousMap.hpp"
#include "HDF5DataSet.hpp"
#include "autodiff/reverse/var/var.hpp"

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
  using Alpha = autodiff::var;
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
  // encapsulates data required to sample a value of beta using proper
  // orthogonal decomposition; limited to a given incident energy range
  struct BetaPartition {
    // constructs a BetaPartition from a `partition` node
    BetaPartition(const pugi::xml_node& partition_node);
    // evaluates beta using proper orthogonal decomposition and linear
    // interpolation in T
    Beta Evaluate(const size_t cdf_index, const size_t E_index, Temperature T)
        const noexcept;
    // contains CDF modes
    const HDF5DataSet<2> CDF_modes;
    // contains singular values
    const HDF5DataSet<1> singular_values;
    // contains energy and temperature modes
    const HDF5DataSet<3> E_T_modes;
  };
  // encapsulates data required to sample a value of alpha using proper
  // orthogonal decomposition; limited to a given beta range
  struct AlphaPartition {
    // constructs an AlphaPartition from a `partition` node
    AlphaPartition(const pugi::xml_node& partition_node);
    // evaluates alpha using proper orthogonal decomposition
    Alpha Evaluate(
        const size_t cdf_index, const size_t local_beta_index,
        const size_t T_index) const;
    // samples a value of alpha using quadratic interpolation in CDF then
    // linear interpolation in T
    Alpha Sample(
        Particle& p, const size_t b_i, const Beta b, const Temperature T,
        const Alpha a_cutoff, const Real awr) const noexcept;
    // contains CDF modes
    const HDF5DataSet<2, true> CDF_modes;
    // contains singular values
    const HDF5DataSet<1, true> singular_values;
    // contains beta and temperature modes
    const HDF5DataSet<3, true> beta_T_modes;
    // total number of elements in this partition (C++ Core Guidelines C.131)
    const size_t size;
  };
  // identifies optimal set of PDFs for a piecewise quadratic fit
  static std::vector<autodiff::var> GetDerivatives(
      const std::vector<autodiff::var>& xs, const std::vector<Real>& ys) noexcept;
  // evaluates a piecewise quadratic function at a point
  static autodiff::var EvaluateQuadratic(
      const std::vector<autodiff::var>& xs, const std::vector<Real>& ys,
      const std::vector<autodiff::var>& fs, Real x) noexcept;
  // performs the inverse of EvaluateQuadratic. Assumes that data is strictly
  // increasing.
  static autodiff::var SolveQuadratic(
      const std::vector<autodiff::var>& xs, const std::vector<Real>& ys,
      const std::vector<autodiff::var>& fs, autodiff::var y) noexcept;
  // evaluates the inelastic scattering cross section.
  MacroscopicCrossSection
  EvaluateInelastic(const size_t E_index, const size_t T_index) const noexcept;
  // Samples an outgoing energy. Requires Particle energy is strictly below
  // ThermalScattering::cutoff_energy. Uses histogram interpolation in PDF.
  Beta SampleBeta(Particle& p, ContinuousEnergy E, Temperature T) const;
  // samples an outgoing cosine given an outgoing energy
  Real SampleAlpha(
      Particle& p, const Beta& b, ContinuousEnergy E, Temperature T) const;
  // majorant cross section
  const ContinuousMap<ContinuousEnergy, MicroscopicCrossSection> majorant;
  // total scattering cross section POD coefficients
  const HDF5DataSet<2> scatter_xs_T;
  const HDF5DataSet<1> scatter_xs_S;
  const HDF5DataSet<2> scatter_xs_E;
  // beta proper orthogonal decomposition coefficients
  const std::vector<BetaPartition> beta_partitions;
  // concatenated vector of all incident energies found in each partition
  const std::vector<ContinuousEnergy> Es;
  // A vector of one-past-the-last energy index for each partition
  const std::vector<size_t> beta_partition_E_ends;
  // alpha proper orthogonal decomposition coefficients
  const std::vector<AlphaPartition> alpha_partitions;
  // concatenated vector of all betas found in each partition_ends
  const std::vector<Beta> betas;
  // identifies the global beta index where each partition begins
  const std::vector<size_t> alpha_partition_beta_begins;
  // start index of each alpha partition for sensitivity analysis
  const std::vector<size_t> alpha_partition_offsets;
  // maximum value of beta which can be sampled
  const Beta beta_cutoff;
  // maximum value of alpha which can be sampled
  const Alpha alpha_cutoff;
  // the target Nuclide
  const Nuclide& target;
  // neutrons above the cutoff energy do not get thermal scattering treatment
  const ContinuousEnergy cutoff_energy = scatter_xs_E.axes.at(0).back();
};
