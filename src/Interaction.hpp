#pragma once

#include "BasicTypes.hpp"
#include "Particle.hpp"

#include <map>
#include <memory>

namespace pugi {
class xml_node;
}

/// @brief Models the interaction between a Particle and a Nuclide
/// @details The polymorphism here shall be where multigroup and continuous
///          energy cross sections are resolved.
class Interaction {
public:
  /// @brief Associates a Particle::Type with a polymorphic pointer to
  ///        Interaction
  using Map = std::map<Particle::Type, std::unique_ptr<const Interaction>>;
  /// @brief Factory method to create multigroup or continuous Interaction
  ///        from a nuclide node
  /// @returns A `std::map` with Particle::Type as keys and `std::unique_ptr`s
  ///          to the constructed Interaction (C++ Core Guidelines R.30)
  /// @exception std::runtime_error Particle declared in `general/particles`
  ///            node not found in `nuclide` node.
  static Map Create(const pugi::xml_node& nuclide_node);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Interaction() noexcept;
  /// @brief Returns the majorant cross section for a given Particle
  virtual MicroscopicCrossSection
  GetMajorant(const Particle& p) const noexcept = 0;
  /// @brief Returns the total cross section for a given Particle
  virtual MicroscopicCrossSection
  GetTotal(const Particle& p) const noexcept = 0;
  /// @brief Interact with a Particle, updating its state
  virtual void Interact(Particle& p) const noexcept = 0;
};
