#pragma once

#include "BasicTypes.hpp"
#include "ContinuousEvaluation.hpp"
#include "ContinuousMap.hpp"

#include <memory>
#include <optional>
#include <vector>

namespace Estimator {
class Proxy;
}
namespace pugi {
class xml_node;
}
class Nuclide;
class Particle;
class ThermalScattering;

/// @brief Abstract interface for reactions which update the state of a
///        Particle
class ContinuousReaction {
public:
  /// @brief Factory method to create a new ContinuousReaction from an XML
  ///        document and target Nuclide
  static std::unique_ptr<const ContinuousReaction> Create(
      const pugi::xml_node& reaction_node, const Nuclide& target,
      const std::optional<ThermalScattering>& tnsl);
  /// @brief Constructs ContinuousReaction from a reaction node of an XML
  ///        document
  ContinuousReaction(
      const pugi::xml_node& reaction_node, const Nuclide& target);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~ContinuousReaction() noexcept;
  /// @brief Returns a cross section that is at least greater than any cross
  ///        section that will be encountered in the current Cell
  /// @details Cross sections are currently given at a single, user-specified
  ///          temperature.
  /// @todo Add warning the first time the requested temperature does not match
  ///       evaluated temperature
  virtual MicroscopicCrossSection
  GetCellMajorant(const Particle& p) const noexcept;
  /// @brief Interact with a Particle, updating its state
  virtual void Interact(
      Particle& p,
      std::vector<Estimator::Proxy>& estimator_proxies) const noexcept = 0;

protected:
  /// @brief Reference to target Nuclide
  const Nuclide& target;
  /// @brief Cross section data associated with reaction
  const ContinuousEvaluation evaluation;
};

/// @brief Contains data required to perform a capture interaction
class ContinuousCapture : public ContinuousReaction {
public:
  /// @brief Constructs a ContinuousCapture from a `capture` node of an XML
  ///        document and target Nuclide
  ContinuousCapture(const pugi::xml_node& capture_node, const Nuclide& target);
  /// @brief Captures the Particle, ending its history
  void Interact(Particle& p, std::vector<Estimator::Proxy>& estimator_proxies)
      const noexcept override;
};

/// @brief Contains data required to perform a scatter interaction
/// @details Additional data on scattering methods is available on @ref
///          stream_delegates
class ContinuousScatter : public ContinuousReaction {
public:
  /// @brief Constructs ContinuousScatter from a `scatter` node of an XML
  ///        document and target Nuclide
  /// @details If present, thermal neutron scattering law data will be added
  ContinuousScatter(
      const pugi::xml_node& scatter_node, const Nuclide& target,
      const std::optional<ThermalScattering>& tnsl);
  /// @brief Returns a cross section that is at least greater than any cross
  ///        section that will be encountered in the current Cell
  /// @details Performs temperature adjustments from free gas or thermal
  ///          scattering, if applicable.
  MicroscopicCrossSection
  GetCellMajorant(const Particle& p) const noexcept override;
  /// @brief Scatter the Particle
  void Interact(Particle& p, std::vector<Estimator::Proxy>& estimator_proxies)
      const noexcept override;

private:
  // Returns adjusted free gas scattering cross section for given temperature
  MicroscopicCrossSection
  GetFreeGasScatterAdjustment(const Particle& p, Temperature T) const noexcept;
  // perform free gas scattering
  void ScatterFreeGas(Particle& p) const noexcept;
  // Thermal neutron scattering law S(a,b,T) of InteractionDelegate of Nuclide
  const std::optional<ThermalScattering>& tnsl;
};

/// @brief Contains data required to perform a fission interaction
class ContinuousFission : public ContinuousReaction {
public:
  /// @brief Constructs ContinuousFission from a `fission` node of an XML
  ///        document and target Nuclide
  ContinuousFission(const pugi::xml_node& fission_node, const Nuclide& target);
  /// @brief Induces a fission event, possibly producing secondary particles
  void Interact(Particle& p, std::vector<Estimator::Proxy>& estimator_proxies)
      const noexcept override;

private:
  // Average number of secondary particles produced per fission
  const std::optional<ContinuousMap<ContinuousEnergy, Real>> nubar;
};
