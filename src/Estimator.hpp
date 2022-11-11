#pragma once

#include "BasicTypes.hpp"
#include "Scorable.hpp"
#include "Sensitivity.hpp"

#include <cstddef>
#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <tuple>
#include <vector>

namespace pugi {
class xml_node;
}
class CurrentEstimator;
class World;
class CSGSurface;
class Perturbation;
class Sensitivity;

/// @brief Contains scores, multipliers, and sensitivities for a quantity of
///        interest
class Estimator : public Scorable {
public:
  /// @brief Visitor interface for interacting with an Estimator
  /// @details Classes which desire different behaviors for interacting with
  ///          different Estimator subclasses can derive from this Visitor
  ///          and implement each of the required interfaces. If different
  ///          return types are required, consider templating.
  class Visitor {
  public:
    /// @brief Return type of call to Visit
    using T = std::optional<std::tuple<size_t, Real>>;
    /// @brief Virtual destructor (C++ Core Guidelines C.127)
    virtual ~Visitor() noexcept;
    /// @brief Interface for double dispatch implementations
    /// @returns A tuple to bin index and score to add; nullopt indicates
    ///          nothing to be scored
    virtual T Visit(const CurrentEstimator& e) const noexcept = 0;
  };
  /// @brief Factory method to create new Estimator from an estimator node of
  ///        an XML document
  static std::unique_ptr<Estimator> Create(
      const pugi::xml_node& estimator_node, const World& world,
      const std::vector<std::unique_ptr<const Perturbation>>& perturbations);
  /// @brief Copy constructor. Deep copies Scores and each Sensitivity object.
  Estimator(const Estimator& other) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Estimator() noexcept;
  /// @brief Interface for double dispatch implementations
  virtual Visitor::T Visit(const Visitor& visitor) const noexcept = 0;
  /// @brief Virtual constructor, used for deep copying an Estimator
  virtual std::unique_ptr<Estimator> Clone() const noexcept = 0;
  /// @brief Returns a string suitable for printing
  std::string to_string(const Real total_weight) const noexcept override;
  /// @brief Add scores of other to this including any Sensitivity objects
  Estimator& operator+=(const Estimator& other) noexcept;
  /// @brief Sensitivity objects which depend on this Estimator (C++ Core
  ///        Guidelines C.131)
  const std::vector<std::unique_ptr<Sensitivity>> sensitivities;

protected:
  /// @brief Constructs an Estimator from an XML node
  Estimator(
      const pugi::xml_node& estimator_node,
      const std::vector<std::unique_ptr<const Perturbation>>&
          perturbations) noexcept;
};

/// @brief A lightweight object for accumulating Estimator scores and any
///        associated Sensitivity scores for a single history
/// @details During a call to Driver::Transport(), scores are produced by a
///          single Particle object's random walk. However, we cannot
///          immediately score an Estimator once that Particle has died. The
///          Particle might have produced secondaries and those secondaries
///          must complete their random walks too. Moreover, the score cannot
///          be squared (for estimating uncertainties) until all scoring for
///          the history is finished.
class EstimatorProxy {
public:
  /// @brief Constructs a EstimatorProxy for an existing Estimator
  EstimatorProxy(Estimator& original) noexcept;
  /// @brief Calls Estimator::Visit for scoring
  void Visit(const Estimator::Visitor& visitor) noexcept;
  /// @brief Adds the pending scores to the actual Estimator
  void CommitHistory() const noexcept;

private:
  // Dependent SensitivtyProxy objects
  std::vector<SensitivityProxy> sensitivity_proxies;
  // map from indices in flattened ParticleBin index to pending scores
  std::map<size_t, Real> pending_scores;
  // Reference to original Estimator
  Estimator& original;
};

/// @brief Scores when a Particle crosses a given CSGSurface
class CurrentEstimator : public Estimator {
public:
  /// @brief Constructs a CurrentEstimator from a `current` node
  CurrentEstimator(
      const pugi::xml_node& current_estimator_node, const World& world,
      const std::vector<std::unique_ptr<const Perturbation>>& perturbations);
  /// @brief Copy constructor; deep copies each sensitivity Estimator in other
  CurrentEstimator(const CurrentEstimator& other) noexcept;
  /// @brief Interface for double dispatch implementations
  Visitor::T Visit(const Visitor& visitor) const noexcept final;
  /// @brief Returns a pointer to a new instance of CurrentEstimator
  std::unique_ptr<Estimator> Clone() const noexcept override;
  /// @brief CSGSurface which this estimator is associated with (C++ Core
  ///        Guidelines C.131). Pointer makes comparison fast.
  const std::shared_ptr<const CSGSurface> surface;
};

/// @brief Collects multiple Estimator objects
class EstimatorSet {
public:
  /// @brief Constructs an EstimatorSet from an XML document
  /// @param estimators_node `estimators` node of an XML document
  /// @param world Used to reference CSGSurface or Cell objects for estimators
  /// @param perturbations Used to associate Estimator objects which have a
  ///        `sensitivities` node with a Perturbation
  /// @param total_weight Used for normalizing Estimator scores
  EstimatorSet(
      const pugi::xml_node& estimators_node, const World& world,
      const std::vector<std::unique_ptr<const Perturbation>>& perturbations,
      const Real total_weight);
  /// @brief Copy constructor; deep copies each Estimator in other
  EstimatorSet(const EstimatorSet& other) noexcept;
  /// @brief Return a collection of EstimatorSet::Proxy
  std::vector<EstimatorProxy> CreateProxies() const noexcept;
  /// @brief Return a reference to Estimator by name
  /// @todo Template along with PerturbationSet::FindPerturbationByName
  const Estimator& FindEstimatorByName(const std::string& name) const;
  /// @brief Returns a string suitable for printing
  std::string to_string() const noexcept;
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
