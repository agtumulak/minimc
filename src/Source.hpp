#pragma once

#include "BasicTypes.hpp"
#include "Particle.hpp"
#include "Point.hpp"

#include <memory>
#include <variant>

namespace pugi {
class xml_node;
}

/// @brief Template for distributions which can be sampled with an RNG
/// @tparam T Returned type of the distribution
template <typename T> class Distribution {
public:
  /// @brief Factory method to create a new Distribution from an XML document
  /// @returns A `std::unique_ptr` to the constructed Distribution (C++ Core
  ///          Guidelines R.30)
  /// @param property_node Either a `position`, `direction`, `energy`, or
  ///        `particletype` node
  static std::unique_ptr<const Distribution<T>>
  Create(const pugi::xml_node& property_node);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Distribution() noexcept;
  /// @brief Samples the distribution
  virtual T Sample(RNG& rng) const noexcept = 0;
};

/// @brief Distribution that always returns the same value
template <typename T> class ConstantDistribution : public Distribution<T> {
public:
  /// @brief Constructs the distribution with a constant value
  ConstantDistribution(const T& constant);
  /// @brief Returns the constant value
  T Sample(RNG& rng) const noexcept override;

private:
  const T constant;
};

/// @brief Returns a point isotropically sampled on the unit sphere
class IsotropicDistribution : public Distribution<Direction>{
  /// @brief Returns an isotropically sampled Direction
  Direction Sample(RNG& rng) const noexcept override;
};

/// @brief Returns an isotropic flux with respect to a given Direction
/// @details An isotropic flux is distributed as @f$ p_{\mu}(\mu) = 2\mu @f$
///          in @f$ [0, 1) @f$ where @f$ \mu @f$ is the cosine of the angle
///          between the sampled Direction and the reference Direction.
class IsotropicFlux : public Distribution<Direction> {
public:
  /// @brief Constructs the distribution with a reference Direction
  IsotropicFlux(const pugi::xml_node& isotropic_flux_node) noexcept;
  /// @brief Returns an isotropic flux
  Direction Sample(RNG& rng) const noexcept override;

private:
  const Direction reference;
};

/// @brief Models a random source of Particle objects
class Source {
public:
  /// @brief Constructs a Source from a `fixedsource` or `initialsource` node
  ///        of an XML document
  /// @details `SourceDistributionType` is a `complexType` defined in the
  ///          minimc XML schema
  Source(const pugi::xml_node& source_node);
  /// @brief Samples a source Particle
  Particle Sample(RNG::result_type seed) const noexcept;

private:
  std::unique_ptr<const Distribution<Point>> position;
  std::unique_ptr<const Distribution<Direction>> direction;
  std::unique_ptr<const Distribution<Energy>> energy;
  std::unique_ptr<const Distribution<Particle::Type>> particle_type;
};
