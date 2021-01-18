#pragma once

#include "BasicTypes.hpp"
#include "Estimator.hpp"
#include "Particle.hpp"
#include "World.hpp"

#include <vector>

/// @brief Represents an independent realization of a history, given an initial
///        Particle
class History {
public:
  /// @brief Returns statistics about a single history
  struct Outcome {
    /// @brief Statistics such as collisions, scatter, etc...
    Estimator estimator;
    /// @brief For k-eigenvalue calculations, this returns particles that are
    ///        added to the fisison bank
    std::vector<Particle> banked;
  };
  /// @brief Constructs a History for a source Particle
  History(Particle& source, const World& world);
  /// @brief Sample an entire History and return outcome as an Estimator
  Outcome Transport();

private:
  Particle p;
  RNG rng{p.seed};
  const World& w;
};
