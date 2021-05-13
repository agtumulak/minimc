#pragma once

#include "BasicTypes.hpp"
#include "NuclearData.hpp"
#include "pugixml.hpp"

#include <random>
#include <string>
#include <vector>

class Particle;

/// @brief Aggregates cross sections for all reactions and related nuclear data
class Nuclide {
public:
  /// @brief Constructs a Nuclide from an XML document
  /// @param root Root node of existing XML document
  /// @param nuclide_name Value of `name` attribute of `nuclide` node in XML
  ///        document
  /// @exception std::runtime_error `nuclide` node with matching `name`
  ///            attribute not found, or incorrect number of entries
  Nuclide(const pugi::xml_node& root, const std::string& nuclide_name);
  /// @brief Returns the total cross section for a given Particle
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept;
  /// @brief Returns the cross section for a given Particle and Reaction
  MicroscopicCrossSection
  GetReaction(const Particle& p, const Reaction r) const noexcept;
  /// @brief Returns the average fission neutron yield for a given Particle
  Real GetNuBar(const Particle& p) const noexcept;
  /// @brief Scatters the Particle and updates its state
  /// @exception std::runtime_error Sampling outgoing Energy failed
  void Scatter(RNG& rng, Particle& p) const;
  /// @brief Fissions the Nuclide and produces secondaries
  std::vector<Particle> Fission(RNG& rng, Particle& p) const noexcept;
  /// @brief Samples a reaction
  Reaction SampleReaction(RNG& rng, const Particle& p) const noexcept;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;

private:
  // Aggregates (polymorphic) NuclearData objects for each Particle::Type
  const NuclearData::Map xs;
};
