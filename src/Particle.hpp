#pragma once

#include "Point.hpp"

/// @brief Particles get transported around
/// @details Particles are described by their position, direction, energy, and
///          an alive flag. Position and direction are represented as Point
///          objects.
/// @note All users of this class should assume that the direction member is
///       normalized to one to avoid extra computation. Maintainers of this
///       class must take care care to keep direction normalized.
class Particle {
public:
  /// @brief Default constructor. Returns a particle at the origin moving in
  ///        the x direction in Group 1.
  Particle() noexcept;
  /// @brief Checks if `*this` is still alive
  explicit operator bool() const noexcept;

private:
  // Position may be anywhere in @f$ \mathbb{R}^3 @f$
  Point position;
  // Direction must be constrained to @f$ \lVert v \rVert = 1 @f$
  Point direction;
  // Energy in continuous energy calculation or group in multigroup calculation
  Energy energy;
  // Flag for determining if this Particle is still alive
  bool alive{true};
};
