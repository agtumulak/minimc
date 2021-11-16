#pragma once

#include "BasicTypes.hpp"
#include "Interaction.hpp"

#include <iosfwd>
#include <string>
#include <vector>

namespace pugi {
class xml_node;
}
class History;
struct State;

/// @brief Aggregates cross sections for all reactions and related nuclear data
class Nuclide {
public:
  /// @brief Constructs a Nuclide from a `nuclide` node
  Nuclide(const pugi::xml_node& nuclide_node);
  /// @brief Returns the majorant cross section for a given State
  MicroscopicCrossSection GetMajorant(const State& s) const noexcept;
  /// @brief Returns the total cross section for a given State
  MicroscopicCrossSection GetTotal(const State& s) const noexcept;
  /// @brief Interact with a particle, updating its State
  void Interact(State& s) const noexcept;
  /// @brief Append secondaries produced by a State onto a History bank
  void AppendSecondaries(
      std::vector<History> bank, const State& state) const noexcept;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;

private:
  // Aggregates (polymorphic) Interaction objects for each Particle::Type
  const Interaction::Map xs;
};
