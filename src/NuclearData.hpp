#pragma once

#include "BasicTypes.hpp"
#include "Particle.hpp"
#include "Reaction.hpp"
#include "pugixml.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

/// @brief Nuclear cross sections in multigroup or continuous energy for a
///        given Nuclide
/// @details The polymorphism here shall be where multigroup and continuous
///          energy cross sections are resolved.
class NuclearData {
public:
  /// @brief Associates a Particle::Type with a polymorphic pointer to
  ///        NuclearData
  using Map = std::map<Particle::Type, std::unique_ptr<const NuclearData>>;
  /// @brief Factory method to create multigroup or continuous NuclearData
  ///        from a nuclide node
  static Map
  Create(const pugi::xml_node& root, const std::string& nuclide_name);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~NuclearData() noexcept;
  /// @brief Returns the total cross section for a given Particle
  virtual MicroscopicCrossSection
  GetTotal(const Particle& p) const noexcept = 0;
  /// @brief Returns the cross section for a given Particle and Reaction
  virtual MicroscopicCrossSection
  GetReaction(const Particle& p, const Reaction r) const noexcept = 0;
  /// @brief Returns the average fission neutron yield for a given Particle
  virtual Real GetNuBar(const Particle& p) const noexcept = 0;
  /// @brief Scatters the Particle and updates its energy and direction
  /// @exception std::runtime_error Most likely cause is that the scattering
  ///            cross section is zero so an outgoing Energy cannot be sampled
  virtual void Scatter(RNG& rng, Particle& p) const = 0;
  /// @brief Fissions the Nuclide and produces secondaries
  virtual std::vector<Particle>
  Fission(RNG& rng, Particle& p) const noexcept = 0;
};
