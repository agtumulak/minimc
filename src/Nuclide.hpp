#pragma once

#include "BasicTypes.hpp"
#include "InteractionDelegate.hpp"
#include "Particle.hpp"

#include <iosfwd>
#include <map>
#include <memory>
#include <string>

namespace pugi {
class xml_node;
}
class EstimatorProxy;

/// @brief Aggregates cross sections for all reactions and related nuclear data
class Nuclide {
public:
  /// @brief Constructs a Nuclide from a `nuclide` node
  Nuclide(const pugi::xml_node& nuclide_node);
  /// @brief Returns the majorant cross section for a given Particle
  MicroscopicCrossSection GetMajorant(const Particle& p) const noexcept;
  /// @brief Returns the total cross section for a given Particle
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept;
  /// @brief Interact with a Particle, updating its state
  void Interact(Particle& p, std::vector<EstimatorProxy>& estimator_proxies)
      const noexcept;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;

private:
  // Aggregates (polymorphic) InteractionDelegate objects for each
  // Particle::Type
  const std::map<Particle::Type, std::unique_ptr<const InteractionDelegate>> xs;
};
