#pragma once

#include "BasicTypes.hpp"
#include "Estimator.hpp"
#include "World.hpp"
#include "pugixml.hpp"

#include <memory>

/// @brief A Driver owns all the data needed to perform radiation transport
class Driver {
public:
  /// @brief Creates objects necessary for a radiation transport
  Driver(const pugi::xml_node& root);
  /// @brief Solves the problem
  virtual void Solve() = 0;

protected:
  /// @brief Accumulates all counts produced by the simulation
  Estimator estimators{};
  /// @brief Global, read-only description of geometric and material properties
  const World world;
  /// @brief Total histories for fixed-source; per-cycle for k-eigenvalue
  const RNG::result_type batchsize;
  /// @brief Number of threads dedicated to particle transport
  const size_t threads;
  /// @brief Used as a parameter for thread-safe classes such as ChunkGiver
  const size_t chunksize;
};
