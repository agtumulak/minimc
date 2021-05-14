#pragma once

#include "BasicTypes.hpp"
#include "Particle.hpp"
#include "Reaction.hpp"
#include "pugixml.hpp"

#include <map>
#include <memory>
#include <string>

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
  static Map
  Create(const pugi::xml_node& root, const std::string& nuclide_name);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Interaction() noexcept;
  /// @brief Returns the total cross section for a given Particle
  virtual MicroscopicCrossSection
  GetTotal(const Particle& p) const noexcept = 0;
  /// @brief Returns the cross section for a given Particle and Reaction
  virtual MicroscopicCrossSection
  GetReaction(const Particle& p, const Reaction r) const noexcept = 0;
  /// @brief Returns the average fission neutron yield for a given Particle
  virtual Real GetNuBar(const Particle& p) const noexcept = 0;
  /// @brief Interact with a Particle, updating its state
  virtual void Interact(Particle& p) const noexcept = 0;
};
