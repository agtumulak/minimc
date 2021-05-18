#pragma once

#include "BasicTypes.hpp"
#include "Interaction.hpp"
#include "Reaction.hpp"
#include "pugixml.hpp"

#include <filesystem>
#include <map>
#include <optional>

class Particle;

/// @brief Contains cross sections which are indexed by continuous energy values
class Continuous : public Interaction {
public:
  /// @brief Constructs continuous energy nuclear data from a particle node of
  ///        an XML document
  Continuous(const pugi::xml_node& particle_node);
  /// @brief Returns the total cross section for a given Particle
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept override;
  /// @brief Returns the cross section for a given Particle and Reaction
  MicroscopicCrossSection
  GetReaction(const Particle& p, const Reaction r) const noexcept override;
  /// @brief Returns the average fission neutron yield for a given Particle
  Real GetNuBar(const Particle& p) const noexcept override;
  /// @brief Interact with a Particle, updating its state
  void Interact(Particle& p) const noexcept override;

private:
  // Continuously maps elements from a domain to a range, provided a limited
  // set of points, by interpolating
  class Map {
  public:
    // Continuous cross sections indexed by one ContinuousEnergy
    using elements_type = std::map<ContinuousEnergy, Real>;
    // Constructs Continuous::Map from a data file
    Map(const std::filesystem::path& datapath);
    // Returns a const reference to the MicroscopicCrossSection at a given
    // ContinuousEnergy
    const Real& at(const ContinuousEnergy e) const noexcept;

  protected:
    elements_type elements;
  };
  // Like Map, but stores elements as the CDF of some random
  // variable. Stores CDF values as keys so that std::map::upper_bound() can be
  // used.
  class CDF : public Map {
  public:
    // Constructs a Continuous::CDF from a data file
    CDF(const std::filesystem::path& datapath);
    // Samples a value from the CDF and returns the sampled key
    elements_type::key_type Sample(RNG& rng) const noexcept;
  };

  using ReactionsMap = std::map<Reaction, Map>;

  // Helper function for reaction cross section construction
  static ReactionsMap CreateReactions(const pugi::xml_node& particle_node);
  // Captures the Particle, killing it
  void Capture(Particle& p) const noexcept;
  // Scatters the Particle and updates its ContinuousEnergy and Direction
  void Scatter(Particle& p) const noexcept;
  /// @brief Fissions the Nuclide and produces secondaries
  void Fission(Particle& p) const noexcept;
  // Average number of secondary particles produced per fission
  const std::optional<Map> nubar;
  // Outgoing energy distribution of fission neutrons
  const std::optional<CDF> chi;
  // Cross section data for each Reaction
  const ReactionsMap reactions;
  // Total cross section provided in nuclear data files
  const Map total;
};
