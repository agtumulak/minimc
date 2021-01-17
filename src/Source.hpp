#pragma once

#include "BasicTypes.hpp"
#include "Particle.hpp"
#include "pugixml.hpp"

#include <random>

class World;

/// @brief Source of particles
class Source {
public:
  /// @brief Constructs a Source from an XML document
  Source(const pugi::xml_node& root);
  /// @brief Samples a Particle.
  /// @brief Returns a Particle with default constructed position, isotropic
  ///        direction, default constructed energy, correctly-set cell,
  ///        default constructed type, and in alive state.
  Particle Sample(std::minstd_rand& rng, const World& w) const noexcept;

private:
  // Helper function to create a set of Particle::Type which this Source will
  // spawn. TODO: Support more than one Particle::Type
  static Particle::Type CreateParticleType(const pugi::xml_node& root) noexcept;
  // Helper function to set default initial energy of spawned Particle
  static Energy CreateDefaultEnergy(const pugi::xml_node& root) noexcept;

  const Particle::Type particle_type;
  const Energy initial_energy;
};
