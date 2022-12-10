#pragma once

#include "BasicTypes.hpp"

#include <map>
#include <memory>
#include <vector>

class Particle;

namespace Perturbation {

namespace IndirectEffect {
class Interface;
}

namespace Sensitivity {

class Interface;

namespace Proxy {

class Visitor;

/// @brief A lightweight object for accumulating sensitivity scores for a
///        single history
class Interface {
public:
  /// @brief Constructs a sensitivity proxy for a sensitivity
  Interface(Sensitivity::Interface& sensitivity) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Interface() noexcept;
  /// @brief Different estimators implement different visitors for interacting
  ///        with a sensitivity proxy
  virtual void Visit(const Visitor& visitor) noexcept = 0;
  /// @brief Each sensitivity proxy must implement a way of adding its
  ///        accumulated scores to the original sensitivity
  virtual void CommitHistory() const noexcept = 0;
  /// @brief Sets this proxy's indirect effects to those of the current particle
  void SetIndirectEffects(const Particle& p) noexcept;
  /// @brief Returns a reference to the indirect effects
  const std::vector<Real>& GetIndirectEfects() const noexcept;

protected:
  /// @brief Care must be taken to reassign indirect effects whenever a new
  ///        Particle is being transported
  std::shared_ptr<const IndirectEffect::Interface> indirect_effect;
  /// @brief Reference to the original sensitivity
  Sensitivity::Interface& sensitivity;
};

/// @brief Accumulates sensitivities for the original sensitivity
class TotalCrossSection : public Interface {
public:
  /// @brief Constructs a sensitivity proxy for a total cross section sensitivty
  TotalCrossSection(Sensitivity::Interface& sensitivity) noexcept;
  /// @brief Implements Interface method
  void Visit(const Visitor& visitor) noexcept final;
  /// @brief Implements Interface method
  void CommitHistory() const noexcept final;
  /// @brief TotalCrossSection only perturbs one parameter so indices correspond
  ///        to bin index
  std::map<BinIndex, Score> pending_scores;
};

}; // namespace Proxy
}; // namespace Sensitivity
}; // namespace Perturbation
