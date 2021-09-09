#pragma once

#include "Driver.hpp"
#include "Estimator.hpp"
#include "Source.hpp"
#include "pugixml.hpp"

#include <atomic>
#include <cstddef>

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
  // fixed source from which new Particle objects can be sampled from
  const Source source;
  // Number of histories completed or initiated by all threads. May exceed
  // batchsize since each thread will call it once before ending.
  std::atomic<size_t> histories_elapsed{0};
};
