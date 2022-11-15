#pragma once

#include "BasicTypes.hpp"
#include "Perturbation.hpp" // in lieu of listing all Perturbation subclasses

#include <memory>

class TotalCrossSectionPerturbationIndirectEffect;

/// @brief Interface for storing indirect effects accumulated by a Particle
/// @details A new history requires a new IndirectEffect for accumulating the
///          total indirect effect caused by a Perturbation
class IndirectEffect {
public:
  /// @brief Visitor interface for interacting with an IndirectEffect
  /// @details Classes which desire different behaviors for interacting with
  ///          different IndirectEffect subclasses can derive from this
  ///          Visitor and implement each of the required interfaces. If
  ///          different return types are required, consider templating.
  class Visitor {
  public:
    /// @brief Virtual destructor (C++ Core Guidelines C.127)
    virtual ~Visitor() noexcept;
    /// @brief Interface for double dispatch implementations
    virtual void Visit(TotalCrossSectionPerturbationIndirectEffect&
                           indirect_effect) const noexcept = 0;
  };
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~IndirectEffect() noexcept;
  /// @brief Interface for double dispatch implementations
  virtual void Visit(const Visitor& visitor) noexcept = 0;
  /// @brief Virtual constructor, used for deep copying an IndirectEffect
  virtual std::unique_ptr<IndirectEffect> Clone() const noexcept = 0;
};

/// @brief IndirectEffect for TotalCrossSectionPerturbation
class TotalCrossSectionPerturbationIndirectEffect : public IndirectEffect {
public:
  /// @brief Constructs a TotalCrossSectionPerturbationIndirectEffect from a
  ///        TotalCrossSectionPerturbation
  TotalCrossSectionPerturbationIndirectEffect(
      const TotalCrossSectionPerturbation& perturbation) noexcept;
  /// @brief Calls IndirectEffect::Visit for updating
  ///        TotalCrossSectionPerturbationIndirectEffect::running_total
  void Visit(const Visitor& visitor) noexcept final;
  /// @brief Returns a pointer to a new instance of
  ///        TotalCrossSectionPerturbationIndirectEffect
  std::unique_ptr<IndirectEffect> Clone() const noexcept final;
  /// @brief The perturbation whose indirect effect is being tracked (C++ Core
  ///        Guidelines C.131)
  const TotalCrossSectionPerturbation& perturbation;
  /// @brief Accumulated indirect effect (C++ Core Guidelines C.131)
  Real running_total;
};
