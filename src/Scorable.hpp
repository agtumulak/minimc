#pragma once

#include "BasicTypes.hpp"

#include <cstddef>
#include <iosfwd>
#include <memory>
#include <vector>

namespace pugi {
class xml_node;
}
class ParticleBins;
class Perturbation;

/// @brief Serves as base class for Estimator and Sensitivity
/// @details Accumulates scores then produces mean and standard deviation
///          estimates using sample statistics
class Scorable {
public:
  /// @brief Constructs a Scorable from a `bins` node of an XML document
  /// @details All scores are initialized to zero
  Scorable(const std::string& name, const pugi::xml_node& bins_node) noexcept;
  /// @brief Constructs an Scorable from an existing estimator and
  ///        perturbation
  /// @details All scores are initialized to zero
  Scorable(
      const Scorable& estimator, const Perturbation& perturbation) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Scorable() noexcept;
  /// @brief Returns a string suitable for printing
  virtual std::string to_string(const Real total_weight) const noexcept = 0;
  /// @brief Adds score to given index
  void AddScore(const size_t index, const Real score) noexcept;
  /// @brief Returns all scores; primarily for unit testing
  Real GetScore(const size_t index, const Real total_weight) const noexcept;
  /// @brief Add scores of other to this
  Scorable& operator+=(const Scorable& other) noexcept;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;
  /// @brief Used to map a Particle to a unique index (C++ Core Guidelines
  ///        C.131)
  const std::shared_ptr<const ParticleBins> bins;

protected:
  /// @brief Returns a string suitable for printing
  /// @details Does not return bins since the derived class will decide if bins
  ///          should be printed. For instance, an Estimator may print the bins
  ///          while the bins are omitted for associated Sensitivity objects.
  std::string GetScoreAsString(const Real total_weight) const noexcept;
  /// @brief Flattened array of scores
  std::vector<Real> scores;
  /// @brief Flattened array of square of scores
  std::vector<Real> square_scores;
};
