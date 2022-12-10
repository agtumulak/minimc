#pragma once

#include "BasicTypes.hpp"
#include "Perturbation/Sensitivity/Proxy/Proxy.hpp"

class Particle;

namespace Perturbation {
namespace Sensitivity {
namespace Proxy {
/// @brief Visitor interface for interacting with a sensitivity proxy
/// @details Estimators which desire different behaviors for interacting with
///          different sensitivity proxies can derive from this Visitor and
///          implement each of the required interfaces.
class Visitor {
public:
  /// @brief Constructs a Visitor for scoring the sensitivity of a score
  Visitor(
      const Particle& p, const BinIndex i, const Score s,
      const bool add_to_existing) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Visitor() noexcept;
  /// @brief Interface for double dispatch implementations
  virtual void Visit(TotalCrossSection& e) const noexcept = 0;
  /// @brief Particle which is scoring
  const Particle& particle;
  /// @brief BinIndex where Particle scored
  const BinIndex index;
  /// @brief The unperturbed estimator Score
  const Score score;
  /// @brief Indicates if the visitor should add to an existing index or
  ///        create a new entry in the sensitivity proxy
  const bool add_to_existing;
};
}; // namespace Proxy
}; // namespace Sensitivity
}; // namespace Perturbation
