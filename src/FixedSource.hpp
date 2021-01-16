#pragma once

#include "BasicTypes.hpp"
#include "Estimator.hpp"
#include "Parallel.hpp"
#include "Source.hpp"
#include "World.hpp"
#include "pugixml.hpp"

/// @brief Drives the history of independent Particle objects
class FixedSource {
public:
  /// @brief Constructs a FixedSource calculation from an XML document
  FixedSource(const pugi::xml_node& root);
  /// @brief Run all histories across using a pool of threads working in chunks
  void PoolSolve();

private:
  Estimator StartWorker();
  Estimator global; // Read-only global estimator. Do NOT write to this.
  const World world;
  const Source source;
  const History histories;
  const size_t threads;
  const size_t chunksize;
  parallel::ChunkGiver history_chunks{histories, chunksize};
};
