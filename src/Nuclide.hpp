#pragma once

#include "BasicTypes.hpp"
#include "InteractionDelegate.hpp"
#include "Particle.hpp"

#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <vector>

namespace Estimator {
class Proxy;
}
namespace pugi {
class xml_node;
}
class ThermalScattering;

/// @brief Aggregates cross sections for all reactions and related nuclear data
class Nuclide {
public:
  /// @brief Constructs a Nuclide from a `nuclide` node
  /// @todo Remove calls to xml_node::root in constructor so that Nuclide can
  ///       be constructed without requiring additional context about XML file.
  ///       The intention is to inline as many of the XML files in the test
  ///       directory within the corresponding unit tests without increasing
  ///       line count too much.
  Nuclide(const pugi::xml_node& nuclide_node);
  /// @brief Returns the majorant cross section for a given Particle
  MicroscopicCrossSection GetMajorant(const Particle& p) const noexcept;
  /// @brief Returns the total cross section for a given Particle
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept;
  /// @brief Interact with a Particle, updating its state
  void Interact(Particle& p, std::vector<Estimator::Proxy>& estimator_proxies)
      const noexcept;
  /// @brief Return reference to optional thermal neutron scattering law data
  /// @details Currently only used for assigning members
  /// @exception std::runtime_error Thermal neutron scattering law data
  ///                               undefined for multigroup physics
  const std::optional<ThermalScattering>& GetTNSL() const;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;
  /// @brief Ratio of nuclide mass to neutron mass
  const Real awr;

private:
  // Aggregates (polymorphic) InteractionDelegate objects for each
  // Particle::Type
  const std::map<Particle::Type, std::unique_ptr<const InteractionDelegate>> xs;
};
