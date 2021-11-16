#pragma once

#include "BasicTypes.hpp"
#include "Particle.hpp"
#include "Point.hpp"

#include <memory>

class Cell;
class CSGSurface;
class Nuclide;
class Source;
class World;

/// @brief Represents a single point in phase space along with members which
///        may speed up transport and/or scoring
/// @details Members are ordered to improve alignment with the exception of
///          rng, which is required to sample the initial state.
struct State {
  /// @brief Mutually exclusive events which can occur
  enum class Event {
    birth,
    scatter,
    capture,
    fission,
    surface_cross,
    leak,
  };
  /// @brief Constructs a State from a Source and RNG
  /// @details World is passed to determine the Cell occuped by the initial
  ///          State
  State(
      const World& world, const Source& source, RNG::result_type seed) noexcept;
  /// @brief Constructs a State representing a particle which has streamed a
  ///        given distance and crossed a CSGSurface to enter a new Cell
  /// @details World is passed to determine the new Cell
  State(
      const State& origin, std::shared_ptr<const CSGSurface> surface,
      const Real distance, const Cell& cell) noexcept;
  /// @brief Constructs a State at the moment it collides with a particular
  ///         Nuclide, before it has changed Direction and Energy
  State(const State& origin, const Real distance) noexcept;
  /// @brief Sample a nuclide to collide with at the current point
  std::shared_ptr<const Nuclide>
  SampleNuclide(const Cell& cell, RNG& rng) noexcept;
  /// @brief Sample a random number uniformly in @f$ [0, 1) @f$
  /// @details This is provided as a convenience function because more
  ///          complicated sampling schemes require multiple uniformly
  ///          distributed random numbers
  Real Sample() noexcept;
  /// @brief Returns true if the State can still transition to a new State
  bool IsAlive() const noexcept;
  /// @brief %State of the random number generator
  RNG rng;
  /// @brief Position @f$ \boldsymbol{x} @f$
  Point position;
  /// @brief Direction of flight @f$ \hat{\boldsymbol{\Omega}} @f$
  Direction direction;
  /// @brief ContinuousEnergy @f$ E @f$ in a continuous energy calculation or
  ///        Group @f$ g @f$ in a multigroup calculation
  Energy energy;
  /// @brief Pointer to surface crossed, if any
  /// @details Not strictly required but this can speed up scoring
  std::shared_ptr<const CSGSurface> surface;
  /// @brief Non-owning pointer to most recent Cell occupied
  /// @details Not strictly required but this can speed up scoring
  const Cell* cell{nullptr};
  /// @brief Pointer to Nuclide collided with, if any
  /// @details Not strictly required but this can speed up banking secondaries
  std::shared_ptr<const Nuclide> nuclide;
  /// @brief Flag @f$ r @f$ describing the event that occured
  Event event;
  /// @brief Type @f$ p @f$ of the particle
  Particle particle;
};
