#pragma once

#include "BasicTypes.hpp"
#include "Estimator.hpp"
#include "Point.hpp"

#include <string>
#include <list>

class Cell;
class Nuclide;
class World;

/// @brief The primary entity performing random walks in a World.
/// @details Particles are characterized by their position, direction, energy,
///          type, and an alive flag. The awkward declaration order of member
///          variables is meant to improve alignment.
/// @note All users of this class should assume that the direction member is
///       normalized to one to avoid extra computation. Maintainers of this
///       class must take care care to keep direction normalized.
class Particle {
public:
  /// @brief The result of a Transport call
  struct TransportOutcome {
    /// @brief Adds the result of another transport result to this result
    TransportOutcome& operator+=(TransportOutcome&& rhs) noexcept;
    /// @brief Estimators scored during transport
    Estimator estimator;
    /// @brief Secondary particles banked during transport
    std::list<Particle> banked;
  };
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
  /// @brief Moves the Particle through each state until it dies.
  TransportOutcome Transport(const World& w) noexcept;
  /// @brief Kill the Particle, stopping the history
  void Kill() noexcept;
  /// @brief Moves the particle along its current direction a given distance
  void Stream(const Real distance) noexcept;
  /// @brief Return the current position of the Particle
  const Point& GetPosition() const noexcept;
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

private:
  // Returns the distance the Particle will travel before colliding
  Real SampleCollisionDistance() noexcept;
  // Sample a Nuclide given that the Particle has collided inside its Cell
  const Nuclide& SampleNuclide() noexcept;
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
  ///          Particle initialized with the same member variables and provided
  ///          with the same World should return the same result after calling
  ///          Transport();
  RNG rng;
  /// @brief Type of the Particle (C++ Core Guidelines C.131)
  const Type type{Type::neutron};

private:
  // Flag for determining if this Particle is still alive
  bool alive{true};
};
