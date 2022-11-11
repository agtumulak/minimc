#pragma once

#include "BasicTypes.hpp"
#include "Perturbation.hpp"
#include "Scorable.hpp"

#include <cstddef>
#include <iosfwd>
#include <map>
#include <memory>
#include <vector>

namespace pugi {
class xml_node;
}
class Estimator;

/// @brief Estimates the sensitivity of an Estimator with respect to a
///        Perturbation
class Sensitivity : public Scorable {
public:
  /// @brief Factory method to create a new Sensitivity from a `perturbation`
  ///        child of a `sensitivities` node of an XML document
  /// @exception std::runtime_error Perturbation with given name not found or
  ///            Sensitivity for Estimator with respect to Perturbation not
  ///            supported
  static std::unique_ptr<Sensitivity> Create(
      const pugi::xml_node& sensitivity_perturbation_node,
      const std::vector<std::unique_ptr<const Perturbation>>& perturbations,
      const Estimator& estimator);
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

/// @brief A lightweight object for accumulating Sensitivity scores for a
///        single history
/// @details Semantics similar to EstimatorProxy
class SensitivityProxy {
public:
  /// @brief Constructs a SensitivtyProxy for an existing Sensitivity
  SensitivityProxy(Sensitivity& original) noexcept;
  /// @brief Adds the pending scores to the actual Sensitivity
  void CommitHistory() const noexcept;

private:
  // map from indices in flattened ParticleBin index to pending scores
  std::map<size_t, Real> pending_scores;
  // Reference to original Sensitivity
  Sensitivity& original;
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
};
