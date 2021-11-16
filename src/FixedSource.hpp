#pragma once

#include "Driver.hpp"
#include "Estimator.hpp"
#include "Source.hpp"

#include <atomic>
#include <cstddef>

namespace pugi {
class xml_node;
}
class History;

/// @brief Creates and executes a fixed-source calculation
class FixedSource : public Driver {
public:
  /// @brief Creates objects necessary for a fixed source calculation
  FixedSource(const pugi::xml_node& root);
  /// @brief Spawn workers to work on chunks of history
  EstimatorSet Solve() override;

private:
  // function executed by a worker on a single thread
  EstimatorSet StartWorker();
  // pushes creates a new History for each secondary produced by a History
  void CreateSecondaries(std::vector<History> bank, const History& h);
  // fixed source from which new Particle objects can be sampled from
  const Source source;
  // Number of histories completed or initiated by all threads. May exceed
  // batchsize since each thread will call it once before ending.
  std::atomic<size_t> histories_elapsed{0};
};
