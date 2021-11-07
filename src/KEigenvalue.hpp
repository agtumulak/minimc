#pragma once

#include "BasicTypes.hpp"
#include "Driver.hpp"
#include "Estimator.hpp"
#include "Particle.hpp"

#include <cstddef>
#include <list>
#include <mutex>
#include <optional>
#include <tuple>

namespace pugi {
class xml_node;
}
class Bank;

/// @brief Creates and executes a k-eigenvalue calculation
class KEigenvalue : public Driver {
public:
  /// @brief Creates objects necessary for a k-eigenvalue calculation
  /// @param root Root node of existing XML document
  KEigenvalue(const pugi::xml_node& root);
  /// @brief Solve sequential fixed-source calculations using a fission bank
  ///        passed between cycles
  /// @todo Fix fission bank sorting. It is not deterministically sampling the
  ///       fission bank to populate the source bank. Consider using the fact
  ///       that each Particle has a `secondaries` member. Consider creating a
  ///       SortableParticle bank derived class with a
  ///       `SortableParticle::parent_history` member.
  /// @todo Reimplement implicit fission estimator
  EstimatorSet Solve() override;

private:
  using Cycle = size_t;
  // Function executed by a worker on a single thread
  std::tuple<Bank, EstimatorSet> StartWorker();
  // Returns the next Particle from source_bank in a thread-safe manner
  std::optional<Particle> NextParticle() noexcept;
  // Allows mutually exclusive access to source_bank
  std::mutex source_bank_mutex;
  // Source Particle objects are selected from a single source bank
  std::list<Particle> source_bank{};
  // Last cycle which is inactive
  const Cycle last_inactive;
  // Last cycle which is active. This is the last cycle overall.
  const Cycle last_active;
  // Previous cycle estimated k-effective. The first cycle guesses k = 1.
  Real k{1};
};
