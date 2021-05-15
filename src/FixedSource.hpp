#pragma once

#include "BasicTypes.hpp"
#include "Driver.hpp"
#include "Estimator.hpp"
#include "Parallel.hpp"
#include "Particle.hpp"
#include "Source.hpp"
#include "pugixml.hpp"

/// @brief Creates and executes a fixed-source calculation
class FixedSource : public Driver {
public:
  /// @brief Creates objects necessary for a fixed source calculation
  FixedSource(const pugi::xml_node& root);
  /// @brief Spawn workers to work on chunks of history
  Estimator Solve() override;
  /// @brief Function executed by a worker on a single thread
  Estimator StartWorker();

private:
  // In a fixed-source calculation a single integer (`history`) uniquely
  // determines the history of a particle.
  Particle Sample(RNG::result_type history) const noexcept;

  ChunkGiver chunk_giver{batchsize, chunksize};
  const Source source;
};
