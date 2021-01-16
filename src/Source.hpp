#pragma once

#include "Particle.hpp"

#include <random>

class World;

/// @brief Source of particles
class Source {
public:
  /// @brief Default constructor. Returns default constructed Particles.
  Source() noexcept;
  /// @brief Samples a Particle.
  /// @brief Returns a Particle with default constructed position, isotropic
  ///        direction, default constructed energy, correctly-set cell,
  ///        default constructed type, and in alive state.
  Particle Sample(std::minstd_rand& rng, const World& w) const noexcept;
};
