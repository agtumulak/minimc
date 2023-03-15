#pragma once

#include "Perturbation/IndirectEffect/IndirectEffect.hpp"

namespace Perturbation {
namespace IndirectEffect {
/// @brief Visitor interface for interacting with an IndirectEffect
/// @details Classes which desire different behaviors for interacting with
///          different IndirectEffect subclasses can derive from this
///          IndirectEffectVisitor and implement each of the required
///          interfaces.
class Visitor {
public:
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Visitor() noexcept;
  /// @brief Interface for interacting with TotalCrossSection
  virtual void Visit(TotalCrossSection& indirect_effect) const noexcept;
  /// @brief Interface for interacting with TNSL
  virtual void Visit(TNSL& indirect_effect) const noexcept;
};
} // namespace IndirectEffect
} // namespace Perturbation
