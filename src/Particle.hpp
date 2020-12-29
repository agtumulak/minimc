#pragma once

#include "BasicTypes.hpp"
#include "Point.hpp"

class Cell;
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
  /// @brief Return the current position of the Particle
  const Point& GetPosition() const noexcept;
  /// @brief Return the current direction of the Particle
  const Point& GetDirection() const noexcept;
  /// @brief Returns the current energy of the Particle
  const Energy& GetEnergy() const noexcept;
  /// @brief Returns a reference to the current Cell the Particle is within
  const Cell& GetCell() const;
  /// @brief Sets the current Cell occupied by the Particle
  void SetCell(const Cell& c) noexcept;

private:
  // Position may be anywhere in @f$ \mathbb{R}^3 @f$
  Point position;
  // Direction must be constrained to @f$ \lVert v \rVert = 1 @f$
  Point direction;
  // Energy in continuous energy calculation or group in multigroup calculation
  Energy energy{Group{1}};
  // Non-owning pointer to current const Cell occupied by the Particle
  const Cell* cell{nullptr};

private:
  // Flag for determining if this Particle is still alive
  bool alive{true};
};
