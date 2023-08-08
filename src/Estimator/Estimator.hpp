#pragma once

#include "BasicTypes.hpp"
#include "Bins.hpp"

#include <memory>
#include <string>
#include <vector>

namespace Perturbation {
namespace Sensitivity {
namespace Proxy {
class Visitor;
}
class Interface;
} // namespace Sensitivity
class Interface;
} // namespace Perturbation
namespace pugi {
class xml_node;
}
class CSGSurface;
class Particle;
class ParticleBins;
class World;

namespace Estimator {

class Visitor;

/// @brief Contains scores, multipliers, and sensitivities to perturbations for
///        a quantity of interest
class Interface {
public:
  /// @brief Factory method to create new Interface from an estimator node of
  ///        an XML document
  static std::unique_ptr<Interface> Create(
      const pugi::xml_node& estimator_node, const World& world,
      const std::vector<std::unique_ptr<const Perturbation::Interface>>&
          perturbations);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Interface() noexcept;
  /// @brief Virtual copy constructor, used for deep copying an estimator
  virtual std::unique_ptr<Interface> Clone() const noexcept = 0;
  /// @brief Interface for double dispatch implementations
  virtual Score Visit(const Visitor& visitor) const noexcept = 0;
  /// @brief Each Estimator must define direct effect for a Perturbation
  virtual std::unique_ptr<const Perturbation::Sensitivity::Proxy::Visitor>
  GetSensitivityProxyVisitor(
      const Particle& p, const BinIndex i, const Score s) const noexcept = 0;
  /// @brief Adds score to given index
  void AddScore(const BinIndex i, const Score s) noexcept;
  /// @brief Returns score at index; primarily for unit testing
  Score GetScore(const BinIndex i, const Real total_weight) const noexcept;
  /// @brief Returns a string suitable for printing
  std::string to_string(const Real total_weight) const noexcept;
  /// @brief Add scores of other to this including any Sensitivity objects
  Interface& operator+=(const Interface& other) noexcept;
  /// Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;
  /// Used to map a Particle to a BinIndex (C++ Core Guidelines C.131)
  const std::shared_ptr<const ParticleBins> bins;
  /// @brief The sensitivity of this estimator with respect to zero or more
  ///        perturbations will be estimated
  const std::vector<std::unique_ptr<Perturbation::Sensitivity::Interface>>
      sensitivities;

protected:
  /// @brief Constructs an Interface from an XML node
  /// @details `EstimatorBase` is a `complexType` defined in the minimc XML
  ///          schema
  Interface(
      const pugi::xml_node& estimator_node,
      const std::vector<std::unique_ptr<const Perturbation::Interface>>&
          perturbations) noexcept;
  /// @brief Copy constructor. Deep copies elements of other.
  Interface(const Interface& other) noexcept;

private:
  // Flattened array; indices correspond to BinIndex
  std::vector<Score> sums = std::vector<Score>(bins->size(), 0.);
  // Flattened array; indices correspond to BinIndex
  std::vector<Score> sum_squares = std::vector<Score>(bins->size(), 0.);
};

/// @brief Scores when a Particle crosses a given CSGSurface
class Current : public Interface {
public:
  /// @brief Constructs a Current from a `current` node
  Current(
      const pugi::xml_node& current_node, const World& world,
      const std::vector<std::unique_ptr<const Perturbation::Interface>>&
          perturbations);
  /// @brief Virtual copy constructor
  std::unique_ptr<Interface> Clone() const noexcept final;
  /// @brief Interface for double dispatch implementations
  Score Visit(const Visitor& visitor) const noexcept final;
  /// @brief Implements Interface method
  std::unique_ptr<const Perturbation::Sensitivity::Proxy::Visitor>
  GetSensitivityProxyVisitor(
      const Particle& p, const BinIndex i, const Score s) const noexcept final;
  /// @brief CSGSurface which this estimator is associated with (C++ Core
  ///        Guidelines C.131). Pointer makes comparison fast.
  const std::shared_ptr<const CSGSurface> surface;
};

}; // namespace Estimator
