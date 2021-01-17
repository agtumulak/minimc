#pragma once

#include "NuclearData.hpp"
#include "Particle.hpp"
#include "pugixml.hpp"

#include <filesystem>
#include <random>
#include <map>

/// @brief Contains cross sections which are indexed by continuous energy values
class Continuous : public NuclearData {
public:
  /// @brief Constructs continuous energy nuclear data from a particle node of
  ///        an XML document
  Continuous(const pugi::xml_node& particle_node);
  /// @brief Returns the total cross section for a given Particle
  CrossSection GetTotal(const Particle& p) const noexcept override;
  /// @brief Scatters the Particle and updates its ContinuousEnergy and
  ///        direction
  void Scatter(std::minstd_rand& rng, Particle& p) const noexcept override;
  /// @brief Samples a Reaction
  Reaction SampleReaction(
      std::minstd_rand& rng, const Particle& p) const noexcept override;

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

  private:
    elements_type elements;
  };

  using ReactionsMap = std::map<NuclearData::Reaction, OneDimensional>;

  // Helper function for reaction cross section construction
  static ReactionsMap CreateReactions(const pugi::xml_node& particle_node);
  // Cross section data for each Reaction
  const ReactionsMap reactions;
  // Total cross section provided in nuclear data files
  const OneDimensional total;
};
