#pragma once

#include "BasicTypes.hpp"

#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace pugi {
class xml_node;
}
class Particle;
class ParticleBins;
class Perturbation;
class ScorableProxy;

/// @brief Serves as base class for Estimator and Sensitivity
/// @details Accumulates scores then produces mean and standard deviation
///          estimates using sample statistics
/// @todo Memoize GetScore since sensitivites call it
class Scorable {
  /// @brief ScorableProxy objects will require access to many private members
  ///        of Scorable
  friend ScorableProxy;

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
  /// @brief Returns what a Particle would score
  virtual Real GetScore(const Particle& p) const noexcept = 0;
  /// @brief Returns a string suitable for printing
  virtual std::string to_string(const Real total_weight) const noexcept = 0;
  /// @brief Returns all scores; primarily for unit testing
  Real GetScore(const size_t index, const Real total_weight) const noexcept;
  /// @brief Add scores of other to this
  Scorable& operator+=(const Scorable& other) noexcept;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;

protected:
  /// @brief Returns a string suitable for printing
  /// @details Does not return bins since the derived class will decide if bins
  ///          should be printed. For instance, an Estimator may print the bins
  ///          while the bins are omitted for associated Sensitivity objects.
  std::string GetScoreAsString(const Real total_weight) const noexcept;
  /// @brief Used to map a Particle to a unique index
  const std::shared_ptr<const ParticleBins> bins;
  /// @brief Flattened array of scores
  std::vector<Real> scores;
  /// @brief Flattened array of square of scores
  std::vector<Real> square_scores;
};

/// @brief A lightweight object for accumulating scores for a single history
/// @details During a call to TransportMethod::Transport(), scores are
///          produced by a single Particle object's random walk. However, we
///          cannot immediately score an Estimator once that Particle has
///          died. The Particle might have produced secondaries and those
///          secondaries must complete their random walks too.
class ScorableProxy {
public:
  /// @brief Construct a ScorableProxy for the given Scorable
  ScorableProxy(Scorable& original) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~ScorableProxy() noexcept;
  /// @brief Add score to pending scores
  virtual void Score(const Particle& p) noexcept;
  /// @brief Adds the pending scores to the actual Estimator
  /// @details Squares of each score are also evaluated for estimating
  ///          uncertainties.
  virtual void CommitHistory() const noexcept;

private:
  // map from indices in flattened vector to pending scores
  std::map<size_t, Real> pending_scores;
  // reference to original Scorable
  Scorable& original;
};
