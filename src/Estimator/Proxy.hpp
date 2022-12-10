#pragma once

#include "BasicTypes.hpp"
#include "Perturbation/Sensitivity/Proxy/Proxy.hpp"

#include <map>
#include <memory>
#include <vector>

class Particle;

namespace Estimator {

class Interface;
class Visitor;

/// @brief A lightweight object for accumulating scores for a single history
/// @details During a call to Driver::Transport(), scores are initially
///          produced by the first Particle object's random walk. However, we
///          cannot immediately score an estimator once that Particle has died.
///          The first Particle might have produced secondaries and those
///          secondaries must complete their random walks too. Moreover, the
///          score cannot be squared (for estimating uncertainties) until all
///          scoring for the history is finished.
class Proxy {
public:
  /// @brief Constructs a estimator proxy for an existing estimator
  Proxy(Interface& estimator) noexcept;
  /// @brief Sets the indirect effects of the current Particle to dependent
  ///        sensitivity proxies
  void SetIndirectEffects(const Particle& p) noexcept;
  /// @brief Calls Estimator::Visit for scoring
  void Visit(const Visitor& visitor) noexcept;
  /// @brief Adds the pending scores to the actual Interface
  void CommitHistory() const noexcept;

private:
  // Reference to original estimator
  Interface& estimator;
  // Dependent sensitivities
  std::vector<std::unique_ptr<Perturbation::Sensitivity::Proxy::Interface>>
      sensitivity_proxies;
  // map from indices in flattened ParticleBin index to pending scores
  std::map<BinIndex, Score> pending_scores;
};
}; // namespace Estimator
