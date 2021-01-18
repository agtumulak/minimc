#pragma once

#include "BasicTypes.hpp"
#include "Estimator.hpp"
#include "Parallel.hpp"
#include "Source.hpp"
#include "World.hpp"
#include "pugixml.hpp"

/// @brief Drives the history of independent Particle objects
/// @details Unlike FixedSourceStandalone, some members of FixedSource are
///          references as they are meant to be used for a single cycle of a
///          broader calculation (such as k-eigenvalue).
class FixedSource {
public:
  /// @brief Constructs a FixedSource calculation with existing data members
  /// @param general_node Top-level general settings node of XML document
  /// @param estimators Global estimators
  /// @param world Geometric and material properites used in transport
  /// @param source Spawns new Particle objects
  FixedSource(
      Estimator& estimators, const World& world, const Source& source,
      History histories, size_t threads, size_t chunksize);
  /// @brief Run all histories across using a pool of threads working in chunks
  void PoolSolve();

private:
  Estimator StartWorker();
  Estimator& global; // Global estimator. Do NOT write to from a worker thread.
  const World& world;
  const Source& source;
  const History histories;
  const size_t threads;
  const size_t chunksize;
  parallel::ChunkGiver history_chunks{histories, chunksize};
};

/// @brief Wraps a FixedSource class. Creates all members from XML document.
class FixedSourceStandalone {
public:
  /// @brief Creates objects necessary for a fixed source calculation
  /// @param Root node of existing XML document
  FixedSourceStandalone(const pugi::xml_node& root);
  /// @brief Creates an instance of FixedSource and computes the solution
  void Solve();

private:
  Estimator global{};
  const World world;
  const Source source;
  const History histories;
  const size_t threads;
  const size_t chunksize;
};
