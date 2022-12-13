#pragma once

#include "Driver.hpp"
#include "Estimator/Estimator.hpp"
#include "Source.hpp"

#include <atomic>
#include <cstddef>
#include <filesystem>
#include <vector>

namespace pugi {
class xml_node;
}

/// @brief Creates and executes a fixed-source calculation
class FixedSource : public Driver {
public:
  /// @brief Creates objects necessary for a fixed source calculation
  /// @param root Root node of existing XML document
  /// @param output_filepath Path to save estimator results.
  FixedSource(
      const pugi::xml_node& root,
      const std::filesystem::path& output_filepath = {});
  /// @brief Spawn workers to work on chunks of history
  const std::vector<std::unique_ptr<Estimator::Interface>>& Solve() const final;

private:
  // function executed by a worker on a single thread
  std::vector<std::unique_ptr<Estimator::Interface>>
  StartWorker(std::atomic<size_t>& histories_elapsed) const;
  // fixed source from which new Particle objects can be sampled from
  const Source source;
};
