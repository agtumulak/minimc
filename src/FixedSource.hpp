#pragma once

#include "Driver.hpp"
#include "Estimator.hpp"
#include "Parallel.hpp"
#include "Source.hpp"
#include "pugixml.hpp"

/// @brief Creates and executes a fixed-source calculation
class FixedSource : public Driver {
public:
  /// @brief Creates objects necessary for a fixed source calculation
  FixedSource(const pugi::xml_node& root);
  /// @brief Spawn workers to work on chunks of history
  Estimator Solve() override;

private:
  // function executed by a worker on a single thread
  Estimator StartWorker();
  // returns a subset of all histores in a thread-safe manner
  ChunkGiver chunk_giver{batchsize, chunksize};
  // fixed source from which new Particle objects can be sampled from
  const Source source;
};
