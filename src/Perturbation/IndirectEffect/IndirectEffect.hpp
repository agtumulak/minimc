#pragma once

#include "BasicTypes.hpp"
#include "Perturbation/Perturbation.hpp"

#include <memory>
#include <vector>

namespace Perturbation {
class Interface;
namespace IndirectEffect {

class Visitor;

/// @brief Interface for storing indirect effects accumulated by a Particle
class Interface {
public:
  /// @brief Constructs an indirect effect with given perturbation
  Interface(const Perturbation::Interface& perturbation) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Interface() noexcept;
  /// @brief Interface for double dispatch implementations
  virtual void Visit(const Visitor& visitor) noexcept = 0;
  /// @brief Virtual constructor, used for copying an indirect effect subclass
  virtual std::unique_ptr<Interface> Clone() const noexcept = 0;
  /// @brief The perturbation whose indirect effect is being tracked (C++ Core
  ///        Guidelines C.131)
  const Perturbation::Interface& perturbation;
  /// @brief Each element corresponds to a single parameter being perturbed
  std::vector<Real> indirect_effects =
      std::vector<Score>(perturbation.n_perturbations, 0.);
};

/// @brief Indirect effect for Perturbation::TotalCrossSection
class TotalCrossSection : public Interface {
public:
  /// @brief Constructs indirect effect for a total cross section perturbation
  ///        of a nuclide
  TotalCrossSection(
      const Perturbation::TotalCrossSection& perturbation) noexcept;
  /// @brief Calls IndirectEffect::Visit for updating
  ///        TotalCrossSectionPerturbationIndirectEffect::running_total
  void Visit(const Visitor& visitor) noexcept final;
  /// @brief Returns a new total cross section perturbation instance
  std::unique_ptr<Interface> Clone() const noexcept final;
  /// @brief The nuclide whose total cross section is being perturbed
  /// @details std::shared_ptr is used so std::map::find can be used on
  ///          Material::afracs
  const std::shared_ptr<const Nuclide> nuclide;
};
}; // namespace IndirectEffect
}; // namespace Perturbation
