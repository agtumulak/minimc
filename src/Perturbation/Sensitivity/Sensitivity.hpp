#pragma once

#include "BasicTypes.hpp"
#include "Bins.hpp"
#include "Perturbation/Perturbation.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace pugi {
class xml_node;
}
namespace Estimator {
class Interface;
}

namespace Perturbation {
namespace Sensitivity {

namespace Proxy {
class Interface;
}

/// @brief Estimates the sensitivity of an estimator with respect to a
///        perturbation.
class Interface {
public:
  /// @brief Factory method to create a new sensitivity from a `sensitivity`
  ///        child of a `perturbations` node of an XML document
  /// @exception std::runtime_error Perturbation with given name not found
  static std::unique_ptr<Interface> Create(
      const pugi::xml_node& perturbations_sensitivity_node,
      const std::vector<std::unique_ptr<const Perturbation::Interface>>&
          perturbations,
      const Estimator::Interface& estimator);
  /// @brief Constructs a sensitivity with given parent estimator and
  ///        perturbation
  Interface(
      const Estimator::Interface& estimator,
      const Perturbation::Interface& perturbation) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Interface() noexcept;
  /// @brief Virtual copy constructor, used for deep copying a sensitivity
  virtual std::unique_ptr<Interface> Clone() const noexcept = 0;
  /// @brief Create a Proxy for this sensitivity
  virtual std::unique_ptr<Proxy::Interface> CreateProxy() noexcept = 0;
  /// @brief Returns a string suitable for printing. Does not print bins since
  ///        the parent estimator will print it.
  virtual std::string to_string(const Real total_weight) const noexcept = 0;
  /// @brief Adds score to given index
  void AddScore(const size_t index, const Score s) noexcept;
  /// @brief Returns score at index; primarily for unit testing
  Score GetScore(const size_t i, const Real total_weight) const noexcept;
  /// @brief Adds sensitivities of other to this
  Interface& operator+=(const Interface& other) noexcept;
  /// @brief Unique name generated from Estimator and Perturbation
  const std::string name;
  /// @brief Perturbation whose effect on an Estimator is being estimated
  const Perturbation::Interface& perturbation;

protected:
  /// @brief Original estimator bin; used to determine scoring vector size
  const ParticleBins& bins;
  /// @brief Flattened array; subclasses determine how to interpret indices
  std::vector<Score> sums =
      std::vector<Score>(bins.size() * perturbation.n_perturbations, 0.);
  /// @brief Flattened array; subclasses determine how to interpret indices
  std::vector<Score> sum_squares =
      std::vector<Score>(bins.size() * perturbation.n_perturbations, 0.);
};

/// @brief Estimates the sensitivity with respect to all parameters of a
///        thermal neutron scattering law dataset
class TNSL : public Interface {
public:
  /// @brief Constructs a thermal neutron scattering law perturbation
  ///        sensitivity from an estimator and perturbation
  TNSL(
      const Estimator::Interface& estimator,
      const Perturbation::Interface& perturbation)
  noexcept;
  /// @brief Virtual copy constructor
  std::unique_ptr<Interface> Clone() const noexcept final;
  /// @brief Create a proxy for this sensitivity
  std::unique_ptr<Proxy::Interface> CreateProxy() noexcept final;
  /// @brief Returns a string of thermal neutron scattering law perturbation
  ///        sensitivites suitable for printing
  std::string to_string(const Real total_weight) const noexcept final;
};

}; // namespace Sensitivity
}; // namespace Perturbation
