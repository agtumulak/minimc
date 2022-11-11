#pragma once

#include "BasicTypes.hpp"

#include <iosfwd>
#include <memory>
#include <string>

namespace pugi {
class xml_node;
}
class Estimator;
class Nuclide;
class Sensitivity;
class TotalCrossSectionPerturbation;
class World;

/// @brief Models a change in a system parameter
class Perturbation {
public:
  /// @brief Visitor interface for interacting with a Perturbation
  /// @details Classes which desire different behaviors for interacting with
  ///          different Perturbation subclasses can derive from this Visitor
  ///          and implement each of the required interfaces. If different
  ///          return types are required, consider templating.
  class Visitor {
  public:
    /// @brief Virtual destructor (C++ Core Guidelines C.127)
    virtual ~Visitor() noexcept;
    /// @brief Interface for double dispatch implementations
    virtual Real
    Visit(const TotalCrossSectionPerturbation& p) const noexcept = 0;
  };
  /// @brief Factory method to create a new Perturbation from a perturbation
  ///        node of an XML document
  static std::unique_ptr<const Perturbation>
  Create(const pugi::xml_node& perturbation_node, const World& world) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Perturbation() noexcept;
  /// @brief Interface for double dispatch implementations
  virtual Real Visit(const Visitor& visitor) const noexcept = 0;
  /// @brief Create a Sensitivity estimator from a `perturbation` child of a
  ///        `sensitivities` node of an XML document
  virtual std::unique_ptr<Sensitivity> CreateSensitivity(
      const pugi::xml_node& perturbation_node,
      const Estimator& estimator) const = 0;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;

protected:
  /// @brief Constructs a Perturbation from a perturbation node
  Perturbation(const pugi::xml_node& perturbation_node) noexcept;
};

/// @brief Models a perturbation in a Nuclide total cross section @f$
///        \Sigma_{\text{perturbed}}(E) = \Sigma_{\text{unperturbed}}(E) +
///        \delta \Sigma @f$
class TotalCrossSectionPerturbation : public Perturbation {
public:
  /// @brief Constructs a TotalCrossSectionPerturbation from a `total` node of
  ///        an XML document
  /// @exception std::runtime_error Perturbed Nuclide name not found in World
  TotalCrossSectionPerturbation(
      const pugi::xml_node& total_node, const World& world);
  /// @brief Interface for double dispatch implementations
  Real Visit(const Visitor& visitor) const noexcept final;
  /// @brief Creates sensitivity estimator for a TotalCrossSectionPerturbation
  ///        if it is supported
  std::unique_ptr<Sensitivity> CreateSensitivity(
      const pugi::xml_node& perturbation_node,
      const Estimator& estimator) const override;
  /// @brief Nuclide whose total cross section is being perturbed
  const std::shared_ptr<const Nuclide> nuclide;
};
