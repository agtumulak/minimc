#pragma once

#include "Particle.hpp"
#include "pugixml.hpp"

#include <map>
#include <random>
#include <string>

/// @brief Nuclear cross sections in multigroup or continuous energy for a
///        given Nuclide
/// @details The polymorphism here shall be where multigroup and continuous
///          energy cross sections are resolved.
class NuclearData {
public:
  /// @brief Cross section values are double
  using CrossSection = double;
  /// @brief All possible (mutually-exclusive, so no total) reactions
  ///        regardless of incident particle type
  enum class Reaction {
    capture,
    scatter,
  };
  /// @brief Associates a Particle::Type with a polymorphic pointer to
  ///        NuclearData
  using Map = std::map<Particle::Type, std::unique_ptr<const NuclearData>>;
  /// @brief Factory method to create multigroup or continuous NuclearData
  ///        from a nuclide node
  static Map
  Create(const pugi::xml_node& root, const std::string& nuclide_name);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~NuclearData() noexcept;
  /// @brief Returns the total cross section for a given Particle
  /// @details The polymorphism here shall be where multigroup and continuous
  ///          energy cross sections are resolved.
  virtual CrossSection GetTotal(const Particle& p) const noexcept = 0;
  /// @brief Scatters the Particle and updates its energy and direction
  /// @exception std::runtime_error Most likely cause is that the scattering
  ///            cross section is zero so an outgoing Energy cannot be sampled
  virtual void Scatter(std::minstd_rand& rng, Particle& p) const = 0;
  /// @brief Samples a reaction
  virtual Reaction
  SampleReaction(std::minstd_rand& rng, const Particle& p) const noexcept = 0;

protected:
  /// @brief Helper function to convert from std::string to Reaction
  static Reaction ToReaction(const std::string& name) noexcept;
};
