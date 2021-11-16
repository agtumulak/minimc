#pragma once

#include "BasicTypes.hpp"
#include "History.hpp"
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
class IsotropicDistribution : public Distribution<Direction> {
  /// @brief Returns an isotropically sampled Direction
  Direction Sample(RNG& rng) const noexcept override;
};

/// @brief Models a random source of particles.
/// @details Each Distribution that makes up this Source is pairwise
///          independent. Future Source objects which account for correlations
///          will probably require sampling methods which accept previously
///          sampled variables as a function parameter:
///          @f[
///            p_{X \mid Y} (x \mid y ) \iff \mathtt{f(x, y = optional)}
///          @f]
class Source {
public:
  /// @brief Constructs a Source from a `fixedsource` or `initialsource` node
  ///        of an XML document
  /// @details `SourceDistributionType` is a `complexType` defined in the
  ///          minimc XML schema
  Source(const pugi::xml_node& source_node);
  /// @brief Samples a position
  Point SamplePosition(RNG& rng) const noexcept;
  /// @brief Samples a direction
  Direction SampleDirection(RNG& rng) const noexcept;
  /// @brief Samples an energy
  Energy SampleEnergy(RNG& rng) const noexcept;
  /// @brief Samples a particle
  Particle SampleParticle(RNG& rng) const noexcept;

private:
  const std::unique_ptr<const Distribution<Point>> position;
  const std::unique_ptr<const Distribution<Direction>> direction;
  const std::unique_ptr<const Distribution<Energy>> energy;
  const std::unique_ptr<const Distribution<Particle>> particle;
};
