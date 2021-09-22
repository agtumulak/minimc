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
class Continuous : public Interaction {
public:
  /// @brief Helper function for constructing CE_XS from a pandas HDF5 file
  static ContinuousMap<ContinuousEnergy, MicroscopicCrossSection>
  ReadPandasHDF5(const std::filesystem::path& datapath);
  /// @brief Constructs continuous energy nuclear data from a particle node of
  ///        an XML document
  Continuous(const pugi::xml_node& particle_node);
  /// @brief Returns the majorant cross section for a given Particle
  /// @details Currently used when there is a continuous dependence on
  ///          temperature, returning the largest cross section across all
  ///          temperatures.
  MicroscopicCrossSection GetMajorant(const Particle& p) const noexcept override;
  /// @brief Returns the total cross section for a given Particle
  /// @details This is not guaranteed to be consistent with the sum of all
  ///          mutually exclusive reactions. The total cross section is meant
  ///          to be a user-provided quantity to speed up calculations.
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept override;
  /// @brief Returns the cross section for a given Particle and Reaction
  MicroscopicCrossSection
  GetReaction(const Particle& p, const Reaction r) const noexcept override;
  /// @brief Returns the average fission neutron yield for a given Particle
  Real GetNuBar(const Particle& p) const noexcept override;
  /// @brief Interact with a Particle, updating its state
  void Interact(Particle& p) const noexcept override;

private:
  // Continuous energy cross sections
  using CE_XS = ContinuousMap<ContinuousEnergy, MicroscopicCrossSection>;
  // Helper function for constructing CE_XS from JANIS Web data file. Returns
  // CE_XS::elements_type because ReadJanisWebCDF wraps this function.
  static CE_XS::elements_type
  ReadJanisWeb(const std::filesystem::path& datapath);
  // Helper function for constructing CDF<ContinuousEnergy> from JANIS Web data
  // file
  static CDF<ContinuousEnergy>::elements_type
  ReadJanisWebCDF(const std::filesystem::path& datapath);
  // Helper function for constructing thermal scattering data S(a,b,T) if found
  static std::optional<ThermalScattering>
  ReadPandasSAB(const pugi::xml_node& tsl_node);
  // Helper function for reaction cross section construction
  static std::map<Reaction, CE_XS>
  CreateReactions(const pugi::xml_node& particle_node);
  // Captures the Particle, killing it
  void Capture(Particle& p) const noexcept;
  // Scatters the Particle and updates its ContinuousEnergy and Direction
  void Scatter(Particle& p) const noexcept;
  /// @brief Fissions the Nuclide and produces secondaries
  void Fission(Particle& p) const noexcept;
  // Average number of secondary particles produced per fission
  const std::optional<ContinuousMap<ContinuousEnergy, Real>> nubar;
  // Outgoing energy distribution of fission neutrons
  const std::optional<CDF<ContinuousEnergy>> chi;
  // Neutron thermal scattering law S(a,b,T)
  const std::optional<ThermalScattering> tsl;
  // Cross section data for each Reaction
  const std::map<Reaction, CE_XS> reactions;
  // Total cross section provided in nuclear data files
  const CE_XS total;
  // Atomic weight ratio of target, yes this is duplicated in tsl::awr
  const Real awr;
};
