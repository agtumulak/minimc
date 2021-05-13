#pragma once

#include "BasicTypes.hpp"
#include "Estimator.hpp"
#include "Point.hpp"
#include "Reaction.hpp"

#include <vector>

class Cell;
class Nuclide;
class World;

/// @brief The primary entity performing random walks in a World.
/// @details Particles are characterized by their position, direction, energy,
///          type, and an alive flag. The awkward declaration order of member
///          variables is meant to improve alignment. On Apple clang version
///          12.0.0 the members are as follows:
///          Member        | Size (bytes)  | Total (bytes)
///          ------------- | ------------- | -------------
///          position      |            24 |            24
///          direction     |            24 |            48
///          energy        |            16 |            64
///          cell*         |             8 |            72
///          type          |             4 |            76
///          seed          |             4 |            80
///          alive         |             1 |            81
///          padding       |             7 |            88
/// @note All users of this class should assume that the direction member is
///       normalized to one to avoid extra computation. Maintainers of this
///       class must take care care to keep direction normalized.
class Particle {
public:
  /// @brief The result of a Transport call
  struct TransportOutcome {
    TransportOutcome& operator+=(const TransportOutcome& rhs) noexcept;
    Estimator estimator;
    std::vector<Particle> banked;
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
      const Type type) noexcept;
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
  void SetDirectionIsotropic(RNG& rng) noexcept;
  /// @brief Returns the current energy of the Particle
  const Energy& GetEnergy() const noexcept;
  /// @brief Updates the current energy of the Particle
  void SetEnergy(const Energy& e) noexcept;
  /// @brief Returns a reference to the current Cell the Particle is within
  const Cell& GetCell() const;
  /// @brief Sets the current Cell occupied by the Particle
  void SetCell(const Cell& c) noexcept;

private:
  // Returns the distance the Particle will travel before colliding
  Real SampleCollisionDistance() noexcept;
  // Sample a Nuclide given that the Particle has collided inside its Cell
  const Nuclide& SampleNuclide() noexcept;
  // Samples a Reaction given that the particle has collided with given Nuclide
  Reaction SampleReaction(const Nuclide& nuclide) noexcept;
  // Random number generator
  RNG rng{0};
  // Position may be anywhere in @f$ \mathbb{R}^3 @f$
  Point position{0, 0, 0};
  // Direction must be constrained to @f$ \lVert v \rVert = 1 @f$
  Direction direction{1, 0, 0};
  // Energy in continuous energy calculation or group in multigroup calculation
  Energy energy{Group{1}};
  // Non-owning pointer to current const Cell occupied by the Particle
  const Cell* cell{nullptr};

public:
  /// @brief Type of the Particle (C++ Core Guidelines C.131)
  const Type type{Type::neutron};
  /// @brief Seed used to initialize rng once this Particle is constructed
  /// @details In monte carlo methods, there is a paradoxical requirement that
  ///          results are completely deterministic. In k-eigenvalue
  ///          calculations, a given fixed-source cycle must yield the same
  ///          result every time the simulation is run. To do this, we note
  ///          that the outcome of sampling a history depends entirely on two
  ///          things: the initial state of the source particle, and the
  ///          initial state of the random number generator (rng) used to
  ///          sample the sequence of events that follow.
  RNG::result_type seed {0};

private:
  // Flag for determining if this Particle is still alive
  bool alive{true};
};
