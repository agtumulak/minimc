#pragma once

#include "BasicTypes.hpp"
#include "Point.hpp"
#include "Sensitivity.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <vector>

namespace pugi {
class xml_node;
}
class Bins;
class World;
class CSGSurface;
class Particle;

/// @brief Scores tallies based on the history of a Particle
class Estimator {
public:
  // index in flattened bin vector
  using FlattenedIndex = size_t;
  /// @brief A lightweight object for accumulating Estimator scores for a
  ///        single history
  /// @details During a call to TransportMethod::Transport(), scores are
  ///          produced by a single Particle object's random walk. However, we
  ///          cannot immediately score an Estimator once that Particle has
  ///          died. The Particle might have produced secondaries and those
  ///          secondaries must complete their random walks too.
  class Proxy {
  public:
    /// @brief Construct a Proxy for the given Estimator
    Proxy(Estimator& init) noexcept;
    /// @brief Add score to pending scores
    void Score(const Particle& p) noexcept;
    /// @brief Adds the pending scores to the actual Estimator
    /// @details Squares of each score are also evaluated for estimating
    ///          uncertainties.
    void CommitHistory() const noexcept;

  private:
    // map from index in flattened bin vector to pending score
    std::map<FlattenedIndex, Real> pending_scores;
    // reference to original Estimator
    Estimator& original;
  };
  /// @brief Factory method to create new Estimator from an estimator node of
  ///        an XML document
  static std::unique_ptr<Estimator>
  Create(const pugi::xml_node& estimator_node, const World& world) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Estimator() noexcept;
  /// @brief Virtual constructor, used for deep copying an EstimatorSet
  virtual std::unique_ptr<Estimator> Clone() const noexcept = 0;
  // interface for getting the direct effect due to a perturbation
  virtual Real
  GetDirectEffect(const Particle& p, const Sensitivity& s) const noexcept = 0;
  /// @brief Returns an output string stream suitable for printing
  std::ostringstream GetPrintable(const Real total_weight) const noexcept;
  /// @brief Add scores of other to this
  Estimator& operator+=(const Estimator& other) noexcept;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;

protected:
  /// @brief Constructs an Estimator from an estimator node
  Estimator(const pugi::xml_node& estimator_node) noexcept;

private:
  // Helper visitor to convert Energy to Real
  struct VisitEnergy {
    Real operator()(const ContinuousEnergy& e) const noexcept { return e; }
    Real operator()(const Group& g) const noexcept { return g; }
  };
  // helper function to construct reference direction if binning over cosines
  static std::optional<Direction>
  CreateDirection(const pugi::xml_node& cosine_node) noexcept;
  // Returns index into flattened vector. The returned index correspond to zero
  // offset from a base index. Zero offset the the unmodified score. Higher
  // offsets correspond to sensitivity scores. The offsets are determined by
  // the order of Sensitivity elements in Estimator::sensitivities.
  FlattenedIndex GetBaseIndex(const Particle& p) const noexcept;
  // interface for scoring
  virtual Real GetScore(const Particle& p) const noexcept = 0;
  // interface for getting multipliers
  virtual std::vector<Real>
  GetMultipliers(const Particle& p) const noexcept = 0;
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
    return strides;
  }
  // reference direction against which direction cosines are computed
  const std::optional<Direction> direction;
  // cosine bins, uses direction member variable
  const std::shared_ptr<const Bins> cosine;
  // energy bins
  const std::shared_ptr<const Bins> energy;
  // sensitivities; the first element is always a NoSensitivity
  const std::vector<std::shared_ptr<const Sensitivity>> sensitivities;
  // contains the size of each element to facilitate constant-time lookup
  const std::array<size_t, 3> strides;
  // flattened array of scores corresponding to different bins and multipliers
  std::vector<Real> scores;
  // flattened array of square of scores corresponding to different bins and
  // multipliers
  std::vector<Real> square_scores;
};

/// @brief Scores when a Particle crosses a given CSGSurface
class CurrentEstimator : public Estimator {
public:
  /// @brief Constructs a CurrentEstimator from a `current` node
  CurrentEstimator(
      const pugi::xml_node& current_estimator_node, const World& world);
  /// @brief Returns a pointer to a new instance of CurrentEstimator
  std::unique_ptr<Estimator> Clone() const noexcept override;
  // interface for getting the direct effect due to a perturbation
  Real GetDirectEffect(
      const Particle& p, const Sensitivity& s) const noexcept override;

private:
  // Implements Estimator method
  Real GetScore(const Particle& p) const noexcept override;
  // Implements Estimator method
  virtual std::vector<Real>
  GetMultipliers(const Particle& p) const noexcept = 0;
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
  /// @brief A lightweight class for accumulating scores produced by a single
  ///        history.
  /// @details Constructing an EstimatorSet at the beginning of each history
  ///          may be expensive since many data st
  class Proxy {
  public:
    /// @brief Construct a Proxy for the given EstimatorSet
    Proxy(const EstimatorSet& init) noexcept;
    /// @brief Score to each Estimator::Proxy
    void Score(const Particle& p) noexcept;
    /// @brief Calls Estimator::Proxy::CommitHistory() on each Estimator::Proxy
    void CommitHistory() const noexcept;

  private:
    // contains all estimator proxies
    std::vector<Estimator::Proxy> estimator_proxies;
  };
  /// @brief Constructs an EstimatorSet from an XML document
  /// @param estimators_node `estimators` node of an XML document
  /// @param world Used to reference CSGSurface or Cell objects for estimators
  /// @param total_weight Used for normalizing Estimator scores
  EstimatorSet(
      const pugi::xml_node& estimators_node, const World& world,
      const Real total_weight);
  /// @brief Copy constructor; deep copies each Estimator in other
  EstimatorSet(const EstimatorSet& other) noexcept;
  /// @brief Return an EstimatorSet::Proxy
  Proxy GetProxy() const noexcept;
  /// @brief Add scores from each Estimator of other to this
  /// @exception An Estimator in this EstimatorSet was not found in the other
  ///            EstimatorSet
  EstimatorSet& operator+=(const EstimatorSet& other);

private:
  // each EstimatorSet contains its own set of pointers to Estimator
  const std::vector<std::unique_ptr<Estimator>> estimators;
  // total weight for normalizing score
  const Real total_weight;
};
