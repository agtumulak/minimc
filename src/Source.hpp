#pragma once

#include "Particle.hpp"

/// @brief Source of particles
class Source {
public:
  /// @brief Default constructor. Returns default constructed Particles.
  Source() noexcept;
  /// @brief Samples a Particle.
  Particle Sample() const noexcept;
};
