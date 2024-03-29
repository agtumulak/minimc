#pragma once

#include "BasicTypes.hpp"
#include "ContinuousEvaluation.hpp"
#include "ContinuousMap.hpp"
#include "ThermalScattering.hpp"

#include <memory>
#include <optional>

namespace pugi {
class xml_node;
}
class Particle;

/// @brief Abstract interface for reactions which update the state of a
///        Particle
class ContinuousReaction {
public:
  /// @brief Factory method to create a new ContinuousReaction from an XML
  ///        document
  static std::unique_ptr<const ContinuousReaction>
  Create(const pugi::xml_node& reaction_node);
  /// @brief Constructs ContinuousReaction from a reaction node of an XML
  ///        document
  ContinuousReaction(const pugi::xml_node& reaction_node);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~ContinuousReaction() noexcept;
  /// @brief Returns true if the reaction's cross section will change the
  ///        total cross section
  /// @details Used when evaluation of the correct total cross section needs
  ///          more information about the reaction. Currently used by
  ///          ContinuousScatter to signal if thermal scattering exists,
  ///          requiring an adjustment of the total cross section, even if the
  ///          total cross section was processed at the correct temperature.
  virtual bool ModifiesTotal(const Particle&) const noexcept;
  /// @brief Returns the largest cross section that may be found within the
  ///        current Cell
  /// @details Currently wraps GetCrossSection
  virtual MicroscopicCrossSection
  GetMajorant(const Particle& p) const noexcept;
  /// @brief Returns temperature-adjusted cross section for the reaction
  /// @details Currently returns the cross section without considering
  ///          temperature
  /// @todo Add warning the first time the requested temperature does not match
  ///       evaluated temperature
  virtual MicroscopicCrossSection
  GetCrossSection(const Particle& p) const noexcept;
  /// @brief Interact with a Particle, updating its state
  virtual void Interact(Particle& p) const noexcept = 0;

protected:
  /// @brief Cross section data associated with reaction
  const ContinuousEvaluation evaluation;
};

/// @brief Contains data required to perform a capture interaction
class ContinuousCapture : public ContinuousReaction {
public:
  /// @brief Constructs a ContinuousCapture from a `capture` node of an XML
  ///        document
  ContinuousCapture(const pugi::xml_node& capture_node);
  /// @brief Captures the Particle, ending its history
  void Interact(Particle& p) const noexcept override;
};

/// @brief Contains data required to perform a scatter interaction
/// @details Additional data on scattering methods is available on @ref
///          transport_methods
class ContinuousScatter : public ContinuousReaction {
public:
  /// @brief Constructs ContinuousScatter from a `scatter` node of an XML
  ///        document, loading thermal scattering data if found
  ContinuousScatter(const pugi::xml_node& scatter_node);
  /// @brief Returns true if thermal scattering is present
  bool ModifiesTotal(const Particle& p) const noexcept override;
  /// @brief Returns largest scattering cross section that may be found within
  ///        the current Cell.
  /// @details Used when there is a continuous dependence on temperature.
  ///          Performs temperature adjustments from free gas or thermal
  ///          scattering, if applicable.
  MicroscopicCrossSection
  GetMajorant(const Particle& p) const noexcept override;
  /// @brief Returns appropriately-adjusted scattering cross section
  /// @details Thermal scattering is assumed to encompass both elastic and
  ///          inelastic scattering. Performs temperature adjustments from free
  ///          gas or thermal scattering, if applicable.
  MicroscopicCrossSection
  GetCrossSection(const Particle& p) const noexcept override;
  /// @brief Scatter the Particle
  /// @todo Compute @f$ \mu @f$ more directly in free gas scattering. Find
  ///       analytic expression to avoid dot product between incident neutron
  ///       velocity and outgoing neutron velocity when calling
  ///       Particle::Scatter.
  void Interact(Particle& p) const noexcept override;

private:
  // Helper function for constructing thermal scattering data S(a,b,T) if found
  static std::optional<ThermalScattering>
  ReadPandasSAB(const pugi::xml_node& ts_node);
  // Returns true if free gas scattering adjustments are applicable. Because
  // downscattering is always possible if the target is hydrogen-1, free gas
  // adjustments are always made if the atomic weight ratio is less than one.
  bool IsFreeGasScatteringValid(
      const Particle& p, const Temperature& T) const noexcept;
  // Returns adjusted free gas scattering cross section for given temperature
  MicroscopicCrossSection
  GetFreeGasScatterAdjustment(const Particle& p, Temperature T) const noexcept;
  // Neutron thermal scattering law S(a,b,T)
  const std::optional<ThermalScattering> tsl;
  // Atomic weight ratio of target, yes this is duplicated in tsl::awr
  const Real awr;
};

/// @brief Contains data required to perform a fission interaction
class ContinuousFission : public ContinuousReaction {
public:
  /// @brief Constructs ContinuousFission from a `fission` node of an XML
  ///        document
  ContinuousFission(const pugi::xml_node& fission_node);
  /// @brief Induces a fission event, possibly producing secondary particles
  void Interact(Particle& p) const noexcept override;

private:
  // Average number of secondary particles produced per fission
  const std::optional<ContinuousMap<ContinuousEnergy, Real>> nubar;
};
