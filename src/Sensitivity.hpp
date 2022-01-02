#pragma once

#include "BasicTypes.hpp"
#include "Perturbation.hpp"
#include "Scorable.hpp"

#include <iosfwd>
#include <memory>

namespace pugi {
class xml_node;
}
class Estimator;
class Particle;

/// @brief Estimates the sensitivity of an Estimator with respect to a
///        Perturbation
class Sensitivity : public Scorable {
public:
  /// @brief Factory method to create a new Sensitivity from a `perturbation`
  ///        child of a `sensitivities` node of an XML document
  /// @exception Sensitivity for Estimator with respect to Perturbation not
  ///            supported
  static std::unique_ptr<Sensitivity> Create(
      const pugi::xml_node& sensitivity_perturbation_node,
      const PerturbationSet& perturbations, const Estimator& estimator);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Sensitivity() noexcept;
  /// @brief Virtual constructor, used for deep copying a Sensitivity
  virtual std::unique_ptr<Sensitivity> Clone() const noexcept = 0;
  /// @brief Returns a string suitable for printing. Does not return bins since
  ///        the parent Estimator will print it.
  std::string to_string(const Real total_weight) const noexcept;
  /// @brief Add scores of other to this
  Sensitivity& operator+=(const Sensitivity& other) noexcept;

protected:
  /// @brief Constructs a Sensitivity with given parent Estimator and
  ///        Perturbation
  Sensitivity(
      const Estimator& estimator, const Perturbation& perturbation) noexcept;
  /// @brief Estimator whose sensitivity to a Perturbation is being estimated
  const Estimator& estimator;
  /// @brief Perturbation whose effect on an Estimator is being estimated
  const Perturbation& perturbation;
};

/// @brief Estimates the sensitivity of a CurrentEstimator with respect to a
///        TotalCrossSectionPerturbation
class CurrentTotalCrossSectionSensitivity : public Sensitivity {
public:
  /// @brief Constructs a CurrentTotalCrossSectionSensitivity from an Estimator
  ///        and Perturbation
  CurrentTotalCrossSectionSensitivity(
      const Estimator& estimator, const Perturbation& perturbation) noexcept;
  /// @brief Returns a pointer to a new instance of
  ///        CurrentTotalCrossSectionSensitivity
  std::unique_ptr<Sensitivity> Clone() const noexcept override;
  /// @brief Implements Scorable method
  Real GetScore(const Particle& p) const noexcept override;
};
