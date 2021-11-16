#pragma once

#include "BasicTypes.hpp"
#include "Point.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace pugi {
class xml_node;
}
class Bins;
class History;
class World;
class CSGSurface;

/// @brief Scores tallies based on the result of a History
/// @todo Add square of scores. This will require <em>the square of the sum of
///       scores produced by one history</em>, <b>not</b> the sum of the square
///       of scores produced by one history.
class Estimator {
  /// @brief Writes contents of Estimator to ostream
  friend std::ostream&
  operator<<(std::ostream& os, const Estimator& e) noexcept;

public:
  /// @brief Factory method to create new Estimator from an estimator node of
  ///        an XML document
  static std::unique_ptr<Estimator>
  Create(const pugi::xml_node& estimator_node, const World& world) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Estimator() noexcept;
  /// @brief Virtual constructor, used for deep copying an EstimatorSet
  virtual std::unique_ptr<Estimator> Clone() const noexcept = 0;
  /// @brief Interface for scoring
  virtual void Score(const History& h) noexcept = 0;
  /// @brief Returns reference to scores vector
  const std::vector<Real>& GetScores() const noexcept;
  /// @brief Normalizes score by total weight
  Estimator& Normalize(Real total_weight) noexcept;
  /// @brief Add scores of other to this
  Estimator& operator+=(const Estimator& other) noexcept;

protected:
  /// @brief Constructs an Estimator from an estimator node
  Estimator(const pugi::xml_node& estimator_node) noexcept;
  /// @brief Return reference to the Bins element where a History would score
  Real& GetScore(const History& h) noexcept;

private:
  // Helper visitor to convert Energy to Real
  struct VisitEnergy {
    Real operator()(const ContinuousEnergy& e) const noexcept { return e; }
    Real operator()(const Group& g) const noexcept { return g; }
  };
  // helper function to construct reference direction if binning over cosines
  static std::optional<Direction>
  CreateDirection(const pugi::xml_node& cosine_node) noexcept;
  // precomputes the stride for each dimension, base case
  template <typename T, typename U>
  static std::array<size_t, 2>
  ComputeStrides(const T&, const U& inner_bin) noexcept {
    return {inner_bin.size(), 1};
  }
  // precompute the strides for each dimension, general case
  template <typename T, typename U, typename... Args>
  static std::array<size_t, sizeof...(Args) + 2> ComputeStrides(
      const T&, const U& middle_bin, const Args&... inner_bins) noexcept {
    std::array<size_t, sizeof...(Args) + 2> strides;
    const auto inner_strides = ComputeStrides(middle_bin, inner_bins...);
    std::copy(
        inner_strides.cbegin(), inner_strides.cend(),
        std::next(strides.begin()));
    strides.front() = middle_bin.size() * inner_strides.front();
  }
  // reference direction against which direction cosines are computed
  const std::optional<Direction> direction;
  // cosine bins, uses direction member variable
  const std::shared_ptr<const Bins> cosine;
  // energy bins
  const std::shared_ptr<const Bins> energy;
  // contains the size of each element to facilitate constant-time lookup
  const std::array<size_t, 2> strides;
  /// @brief Flattened array of scores
  std::vector<Real> scores;
};

/// @brief Scores when a Particle crosses a given CSGSurface
class CurrentEstimator : public Estimator {
public:
  /// @brief Constructs a CurrentEstimator from a `current` node
  CurrentEstimator(
      const pugi::xml_node& current_estimator_node, const World& world);
  /// @brief Returns a pointer to a new instance of CurrentEstimator
  std::unique_ptr<Estimator> Clone() const noexcept override;
  /// @brief Implements Estimator method
  void Score(const History& h) noexcept override;

private:
  // CSGSurface which this estimator is associated with, pointer makes
  // comparison fast
  const std::shared_ptr<const CSGSurface> surface;
};

/// @brief Collects multiple Estimator objects
class EstimatorSet {
  /// @brief Writes contents of Estimator to ostream
  friend std::ostream&
  operator<<(std::ostream& os, const EstimatorSet& e) noexcept;

public:
  /// @brief Copy constructor; deep copies each Estimator in other
  EstimatorSet(const EstimatorSet& other) noexcept;
  /// @brief Constructs an EstimatorSet from an XML document
  /// @param estimators_node `estimators` node of an XML document
  /// @param world Used to reference CSGSurface or Cell objects for estimators
  EstimatorSet(const pugi::xml_node& estimators_node, const World& world);
  /// @brief Return a reference to Estimator by name
  const Estimator& GetEstimator(const std::string& name) const;
  /// @brief Scores each Estimator present in the problem
  void Score(const History& h) noexcept;
  /// @brief Normalizes score of each Estimator by total weight
  EstimatorSet& Normalize(Real total_weight) noexcept;
  /// @brief Add scores from each Estimator of other to this
  EstimatorSet& operator+=(const EstimatorSet& other) noexcept;

private:
  // helper function to deep copy all Estimator objects in other
  static std::map<std::string, std::unique_ptr<Estimator>>
  ConstructAllEstimators(const EstimatorSet& other) noexcept;
  // helper function to construct all Estimator objects from an `estimators`
  // node
  static std::map<std::string, std::unique_ptr<Estimator>>
  ConstructAllEstimators(
      const pugi::xml_node& esitmators_node, const World& world);
  // each EstimatorSet contains its own set of pointers to Estimator
  const std::map<std::string, std::unique_ptr<Estimator>> estimators;
};
