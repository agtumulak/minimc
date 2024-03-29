#pragma once

#include "BasicTypes.hpp"
#include "Interaction.hpp"

#include <iosfwd>
#include <string>

namespace pugi {
class xml_node;
}
class Particle;

/// @brief Aggregates cross sections for all reactions and related nuclear data
class Nuclide {
public:
  /// @brief Constructs a Nuclide from a `nuclide` node
  Nuclide(const pugi::xml_node& nuclide_node);
  /// @brief Returns the majorant cross section for a given Particle
  MicroscopicCrossSection GetMajorant(const Particle& p) const noexcept;
  /// @brief Returns the total cross section for a given Particle
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept;
  /// @brief Interact with a Particle, updating its state
  void Interact(Particle& p) const noexcept;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;

private:
  // Aggregates (polymorphic) Interaction objects for each Particle::Type
  const Interaction::Map xs;
};
