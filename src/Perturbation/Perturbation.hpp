#pragma once

#include <cstddef>
#include <iosfwd>
#include <memory>
#include <string>

namespace Estimator {
class Interface;
}
namespace pugi {
class xml_node;
}
class Nuclide;
class World;

namespace Perturbation {

namespace IndirectEffect {
class Interface;
}
namespace Sensitivity {
class Interface;
}

/// @brief Models a change in a system parameter
/// @details Most perturbations will perturb only a single parameter. However,
///          perturbations <em>can</em> perturb more than one parameter.
///          Encapsulating many perturbations within a Perturbation::Interface
///          is done to avoid repeated evaluations.
class Interface {
public:
  /// @brief Factory method to create a new perturbation from a perturbation
  ///        node of an XML document
  static std::unique_ptr<Interface>
  Create(const pugi::xml_node& perturbation_node, const World& world) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Interface() noexcept;
  /// @brief Each perturbation must implement a Sensitivity::Interface
  virtual std::unique_ptr<Sensitivity::Interface>
  CreateSensitivity(const Estimator::Interface& estimator) const noexcept = 0;
  /// @brief Each perturbation must implement an IndirectEffect::Interface
  virtual std::unique_ptr<IndirectEffect::Interface>
  CreateIndirectEffect() const noexcept = 0;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;
  /// @brief Total number of parameters being perturbed
  const size_t n_perturbations;

protected:
  /// @brief Constructs an Interface from a perturbation node
  /// @details `PerturbationInterface` is a `complexType` defined in the minimc
  ///          XML schema
  Interface(
      const pugi::xml_node& perturbation_node,
      const size_t n_perturbations) noexcept;
};

/// @brief Models a perturbation in a Nuclide microscopic total cross section
///        @f$ \sigma_{\text{perturbed}}(E) = \sigma_{\text{unperturbed}}(E) +
///        \delta \sigma @f$
class TotalCrossSection : public Interface {
public:
  /// @brief Constructs a TotalCrossSection from a `totalxs` node of an XML
  ///        document
  /// @exception std::runtime_error Perturbed Nuclide name not found in World
  TotalCrossSection(const pugi::xml_node& totalxs_node, const World& world);
  /// @brief Returns a Sensitiviity::TotalCrossSection
  std::unique_ptr<Sensitivity::Interface>
  CreateSensitivity(const Estimator::Interface& estimator) const noexcept final;
  /// @brief Returns a IndirectEffect::TotalCrossSection
  std::unique_ptr<IndirectEffect::Interface>
  CreateIndirectEffect() const noexcept final;
  /// @brief Nuclide whose microscopic total cross section is being perturbed
  const std::shared_ptr<const Nuclide> nuclide;
};

}; // namespace Perturbation
