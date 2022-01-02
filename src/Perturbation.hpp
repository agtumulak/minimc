#pragma once

#include "BasicTypes.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace pugi {
class xml_node;
}
class Estimator;
class Nuclide;
class Particle;
class Sensitivity;
class World;

/// @brief Models a change in a system parameter
/// @details A new Perturbation is created for each Particle instance since
///          each Particle will have its own sest of indirect effects which
///          must be accumulated
class Perturbation {
public:
  /// @brief Factory method to create a new Perturbation from a perturbation
  ///        node of an XML document
  static std::unique_ptr<Perturbation>
  Create(const pugi::xml_node& perturbation_node, const World& world) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Perturbation() noexcept;
  /// @brief Create a Sensitivity estimator from a `perturbation` child of a
  ///        `sensitivities` node of an XML document
  virtual std::unique_ptr<Sensitivity> CreateSensitivity(
      const pugi::xml_node& perturbation_node,
      const Estimator& estimator) const = 0;
  /// @brief Returns the indirect effect as a result of performing
  ///        Particle::Stream
  virtual Real
  Stream(const Particle& p, const Real distance) const noexcept = 0;
  /// @brief Returns the indirect effect as a result of performing
  ///        Particle::Scatter
  virtual Real Scatter(
      const Particle& p, const Real& mu, const Energy& e) const noexcept = 0;
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
  /// @brief Creates sensitivity estimator for a TotalCrossSectionPerturbation
  ///        if it is supported
  std::unique_ptr<Sensitivity> CreateSensitivity(
      const pugi::xml_node& perturbation_node,
      const Estimator& estimator) const override;
  /// @brief The total cross section affects the probability of scattering a
  ///        given distance
  Real Stream(const Particle& p, const Real distance) const noexcept override;
  /// @brief The total cross section will not affect the probability of
  ///        scattering with some outgoing angle and energy
  Real
  Scatter(const Particle&, const Real&, const Energy&) const noexcept override;

private:
  // Nuclide whose total cross section is being perturbed
  const std::shared_ptr<const Nuclide> nuclide;
};

/// @brief Aggregates perturbations whose effect on one or more Estimator
///        objects is to be estimated
class PerturbationSet {
public:
  /// @brief Constructs a PerturbationSet from a `perturbations` node of an XML
  ///        document
  PerturbationSet(
      const pugi::xml_node& perturbations_node, const World& world) noexcept;
  /// @brief Returns the Perturbation with the given name
  /// @exception std::runtime_error Perturbation with given name not found
  /// @todo Template along with EstimatorSet::FindEstimatorByName
  const Perturbation& FindPerturbationByName(const std::string& name) const;
  /// @brief Returns a map of indirect effects for use by a Particle
  /// @note Raw pointers are used since a PerturbationSet should outlive any
  ///       Particle instance. Particles are frequently constructed and
  ///       destructed so a std::shared_ptr would update the control block
  ///       quite frequently.
  std::map<const Perturbation*, Real> GetIndirectEffects() const noexcept;

private:
  // PerturbationSet owns each of its Perturbation objects
  const std::vector<std::unique_ptr<const Perturbation>> perturbations;
};
