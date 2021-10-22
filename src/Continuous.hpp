#pragma once

#include "BasicTypes.hpp"
#include "ContinuousMap.hpp"
#include "Interaction.hpp"
#include "Reaction.hpp"
#include "ThermalScattering.hpp"
#include "pugixml.hpp"

#include <filesystem>
#include <map>
#include <optional>

class Particle;

/// @brief Contains cross sections which are indexed by continuous energy values
/// @todo Compute @f$ \mu @f$ more directly in free gas scattering. Find
///       analytic expression to avoid dot product between incident neutron
///       velocity and outgoing neutron velocity when calling Particle::Scatter.
class Continuous : public Interaction {
public:
  /// @brief Constructs continuous energy nuclear data from a particle node of
  ///        an XML document
  Continuous(const pugi::xml_node& particle_node);
  /// @brief Returns the majorant cross section for a given Particle
  /// @details Currently used when there is a continuous dependence on
  ///          temperature, returning the largest cross section across all
  ///          temperatures. Performs temperature adjustments from free gas or
  ///          thermal scattering, if applicable.
  MicroscopicCrossSection GetMajorant(const Particle& p) const noexcept override;
  /// @brief Returns the total cross section for a given Particle
  /// @details This is not guaranteed to be consistent with the sum of all
  ///          mutually exclusive reactions. The total cross section is meant
  ///          to be a user-provided quantity to speed up calculations.
  ///          Performs temperature adjustments from free gas or thermal
  ///          scattering, if applicable.
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept override;
  /// @brief Returns the cross section for a given Particle and Reaction
  MicroscopicCrossSection
  GetReaction(const Particle& p, const Reaction r) const noexcept override;
  /// @brief Returns the average fission neutron yield for a given Particle
  Real GetNuBar(const Particle& p) const noexcept override;
  /// @brief Interact with a Particle, updating its state
  void Interact(Particle& p) const noexcept override;

private:
  // Bundles cross section and temperature data
  struct Evaluation1D {
    // Constructs evaluation from an evaluation node of an XML document
    Evaluation1D(const pugi::xml_node& evaluation_node);
    // Pointwise values of the evaluation
    const ContinuousMap<ContinuousEnergy, MicroscopicCrossSection> xs;
    // The temperature at which the data was evaluated
    const Temperature temperature;
  };
  // Helper function for constructing thermal scattering data S(a,b,T) if found
  static std::optional<ThermalScattering>
  ReadPandasSAB(const pugi::xml_node& tsl_node);
  // Helper function for reaction cross section construction
  static std::map<Reaction, Evaluation1D>
  CreateReactions(const pugi::xml_node& particle_node);
  // Captures the Particle, killing it
  void Capture(Particle& p) const noexcept;
  // Scatters the Particle and updates its ContinuousEnergy and Direction
  void Scatter(Particle& p) const noexcept;
  // Fissions the Nuclide and produces secondaries
  void Fission(Particle& p) const noexcept;
  // Returns true if free gas scattering adjustments are applicable. Because
  // downscattering is always possible if the target is hydrogen-1, free gas
  // adjustments are always made if the atomic weight ratio is less than one.
  bool IsFreeGasScatteringValid(
      const Particle& p, const Temperature& T) const noexcept;
  // Returns adjusted free gas scattering cross section for given temperature
  MicroscopicCrossSection
  GetFreeGasScatterAdjustment(const Particle& p, Temperature T) const noexcept;
  // Average number of secondary particles produced per fission
  const std::optional<ContinuousMap<ContinuousEnergy, Real>> nubar;
  // Neutron thermal scattering law S(a,b,T)
  const std::optional<ThermalScattering> tsl;
  // Cross section data for each Reaction
  const std::map<Reaction, Evaluation1D> reactions;
  // Total cross section provided in nuclear data files
  const Evaluation1D total;
  // Atomic weight ratio of target, yes this is duplicated in tsl::awr
  const Real awr;
};
