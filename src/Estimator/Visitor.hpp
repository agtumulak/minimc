#pragma once

#include "BasicTypes.hpp"
#include "Estimator/Estimator.hpp"

class Particle;

namespace Estimator {

/// @brief Double dispatch visitor interface for interacting with an Estimator
/// @details Classes which desire different behaviors for interacting with
///          different Estimator subclasses can derive from this visitor
///          and implement each of the required interfaces.
class Visitor {
public:
  /// @brief Constructs a visitor using a Particle
  Visitor(const Particle& p) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Visitor() noexcept;
  /// @brief Each subclass must provide an interface for current estimators
  virtual Score Visit(const Current& e) const noexcept = 0;
  /// @brief Reference to the Particle which is scoring; used for accessing
  ///        bin index and indirect effects
  const Particle& particle;
};
}; // namespace Estimator
