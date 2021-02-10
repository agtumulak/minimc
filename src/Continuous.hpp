#pragma once

#include "NuclearData.hpp"
#include "Particle.hpp"
#include "pugixml.hpp"

#include <filesystem>
#include <random>
#include <map>
#include <optional>

/// @brief Contains cross sections which are indexed by continuous energy values
class Continuous : public NuclearData {
public:
  /// @brief Constructs continuous energy nuclear data from a particle node of
  ///        an XML document
  Continuous(const pugi::xml_node& particle_node);
  /// @brief Returns the total cross section for a given Particle
  CrossSection GetTotal(const Particle& p) const noexcept override;
  /// @brief Returns the fission cross section for a given Particle
  CrossSection GetFission(const Particle& p) const noexcept override;
  /// @brief Returns the average fission neutron yield for a given Particle
  Real GetNuBar(const Particle& p) const noexcept override;
  /// @brief Scatters the Particle and updates its ContinuousEnergy and
  ///        direction
  void Scatter(RNG& rng, Particle& p) const noexcept override;
  /// @brief Fissions the Nuclide and produces secondaries
  std::vector<Particle> Fission(RNG& rng, Particle& p) const noexcept override;
  /// @brief Samples a Reaction
  Reaction SampleReaction(RNG& rng, const Particle& p) const noexcept override;

private:
  using elements_type = std::map<ContinuousEnergy, Real>;
  // Continuous cross sections indexed by one ContinuousEnergy
  class OneDimensional {
  public:
    // Constructs Continuous::OneDimensional from a data file
    OneDimensional(const std::filesystem::path& datapath);
    // Returns a const reference to the CrossSection at a given
    // ContinuousEnergy
    const Real& at(const ContinuousEnergy e) const noexcept;

  protected:
    elements_type elements;
  };
  // Like OneDimensional, but stores elements as the CDF of some random
  // variable. Stores CDF values as keys so that std::map::upper_bound() can be
  // used.
  class CDF : public OneDimensional {
  public:
    // Constructs a Continuous::CDF from a data file
    CDF(const std::filesystem::path& datapath);
    // Samples a value from the CDF and returns the sampled key
    elements_type::key_type Sample(RNG& rng) const noexcept;
  };

  using ReactionsMap = std::map<NuclearData::Reaction, OneDimensional>;

  // Helper function for reaction cross section construction
  static ReactionsMap CreateReactions(const pugi::xml_node& particle_node);
  // Average number of secondary particles produced per fission
  const std::optional<OneDimensional> nubar;
  // Outgoing energy distribution of fission neutrons
  const std::optional<CDF> chi;
  // Cross section data for each Reaction
  const ReactionsMap reactions;
  // Total cross section provided in nuclear data files
  const OneDimensional total;
};
