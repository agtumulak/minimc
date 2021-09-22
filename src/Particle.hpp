#pragma once

#include "BasicTypes.hpp"
#include "Point.hpp"

#include <iosfwd>
#include <list>

class Cell;
class Nuclide;

/// @brief The primary entity performing random walks in a World.
/// @details Particles are characterized by their position, direction, energy,
///          type, and an alive flag. The awkward declaration order of member
///          variables is meant to improve alignment.
/// @note Users of this class should assume that Direction is normalized to
///       avoid extra computation.
class Particle {
public:
  /// @brief Affects which cross section data is used during transport, among
  ///        other things
  enum class Type {
    neutron,
    photon,
  };
  /// @brief Helper function to convert from std::string to Type
  static Type ToType(const std::string& name) noexcept;
  /// @brief Member constructor. Explicitly assigns phase-space members.
  Particle(
      const Point& position, const Direction& direction, const Energy& energy,
      const Type type, RNG::result_type seed,
      const Cell* cell = nullptr) noexcept;
  /// @brief Kill the Particle, stopping the history
  void Kill() noexcept;
  /// @brief Moves the particle along its current direction a given distance
  void Stream(const Real distance) noexcept;
  /// @brief Scatters the Particle with an outgoing direction and energy.
  /// @details Scattering is assumed to be azimuthally symmetric.
  /// @param mu The scattering cosine @f$ \mu @f$
  /// @param e The outgoing energy
  void Scatter(const Real& mu, const Energy& e) noexcept;
  /// @brief Return the current position of the Particle
  const Point& GetPosition() const noexcept;
  /// @brief Return the current direction of the Particle
  const Direction& GetDirection() const noexcept;
  /// @brief Sets the direction to a random isotropic direction
  /// @note This should be replaced by a method which accepts scattering cosine
  void SetDirectionIsotropic() noexcept;
  /// @brief Returns the current energy of the Particle
  const Energy& GetEnergy() const noexcept;
  /// @brief Updates the current energy of the Particle
  void SetEnergy(const Energy& e) noexcept;
  /// @brief Returns the Type of the Particle
  Type GetType() const noexcept;
  /// @brief Returns a reference to the current Cell the Particle is within
  const Cell& GetCell() const;
  /// @brief Sets the current Cell occupied by the Particle
  void SetCell(const Cell& c) noexcept;
  /// @brief Banks secondaries produced during transport using an outgoing
  ///        Direction and outgoing Energy
  void Bank(const Direction& direction, const Energy& energy) noexcept;
  /// @brief Sample a random number uniformly in @f$ [0, 1) @f$
  /// @details This is provided as a convenience function because more
  ///          complicated sampling schemes require multiple uniformly
  ///          distributed random numbers
  Real Sample() noexcept;
  /// @brief Sample a Nuclide given that the Particle has collided inside its
  ///        Cell
  const Nuclide& SampleNuclide() noexcept;

private:
  // Secondaries produced
  std::list<Particle> secondaries;
  // Position may be anywhere in @f$ \mathbb{R}^3 @f$
  Point position{0, 0, 0};
  // Direction must be constrained to @f$ \lVert v \rVert = 1 @f$
  Direction direction{1, 0, 0};
  // Energy in continuous energy calculation or group in multigroup calculation
  Energy energy{Group{1}};
  // Non-owning pointer to current const Cell occupied by the Particle
  const Cell* cell{nullptr};

public:
  /// @brief Random number generator (C++ Core Guidelines C.131)
  /// @details This Particle contains its own random number generator. Any
  ///          Particle initialized with the same member variables and passed
  ///          to the same TransportMethod::Transport call should return the
  ///          same TransportMethod::Outcome.
  RNG rng;
  /// @brief Type of the Particle (C++ Core Guidelines C.131)
  const Type type{Type::neutron};
  /// @brief Flag for determining if this Particle is still alive (C++ Core
  ///        Guidelines C.131)
  bool alive{true};
};
