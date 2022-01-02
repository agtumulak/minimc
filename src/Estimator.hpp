#pragma once

#include "BasicTypes.hpp"
#include "Scorable.hpp"

#include <iosfwd>
#include <memory>
#include <vector>

namespace pugi {
class xml_node;
}
class EstimatorProxy;
class EstimatorSetProxy;
class World;
class CSGSurface;
class Particle;
class PerturbationSet;
class Sensitivity;

/// @brief Contains scores, multipliers, and sensitivities for a quantity of
///        interest
class Estimator : public Scorable {
  /// @brief EstimatorProxy objects will require access to many private members
  ///        of Estimator
  friend EstimatorProxy;

public:
  /// @brief Factory method to create new Estimator from an estimator node of
  ///        an XML document
  static std::unique_ptr<Estimator> Create(
      const pugi::xml_node& estimator_node, const World& world,
      const PerturbationSet& perturbations);
  /// @brief Copy constructor. Deep copies Scores and each Sensitivity object.
  Estimator(const Estimator& other) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Estimator() noexcept;
  /// @brief Virtual constructor, used for deep copying an Estimator
  virtual std::unique_ptr<Estimator> Clone() const noexcept = 0;
  /// @brief Returns a string suitable for printing
  std::string to_string(const Real total_weight) const noexcept override;
  /// @brief Add scores of other to this including any Sensitivity objects
  Estimator& operator+=(const Estimator& other) noexcept;

protected:
  /// @brief Constructs an Estimator from an XML node
  Estimator(
      const pugi::xml_node& estimator_node,
      const PerturbationSet& perturbations) noexcept;

private:
  // Sensitivity objects which depend on this Estimator
  std::vector<std::unique_ptr<Sensitivity>> sensitivities;
};

/// @brief A lightweight object for accumulating Estimator scores for a single
///        history
/// @details In addition to scoring its own pending scores, an EstimatorProxy
///          scores pending scores for its Sensitivity objects
class EstimatorProxy : public ScorableProxy {
public:
  /// @brief Construct a Proxy for the given Estimator
  EstimatorProxy(Estimator& original) noexcept;
  /// @brief Add own pending scores as well as pending scores for each
  ///        Sensitivity objects
  void Score(const Particle& p) noexcept override;
  /// @brief Adds the pending scores to the actual Estimator
  void CommitHistory() const noexcept override;

private:
  // Dependent SensitivtyProxy objects; currently do not need anything more
  // derived than ScorableProxy
  std::vector<ScorableProxy> sensitivity_proxies;
};

/// @brief Scores when a Particle crosses a given CSGSurface
class CurrentEstimator : public Estimator {
public:
  /// @brief Constructs a CurrentEstimator from a `current` node
  CurrentEstimator(
      const pugi::xml_node& current_estimator_node, const World& world,
      const PerturbationSet& perturbations);
  /// @brief Copy constructor; deep copies each sensitivity Estimator in other
  CurrentEstimator(const CurrentEstimator& other) noexcept;
  /// @brief Returns a pointer to a new instance of CurrentEstimator
  std::unique_ptr<Estimator> Clone() const noexcept override;
  /// @brief Implements Scorable method
  Real GetScore(const Particle& p) const noexcept override;

private:
  // CSGSurface which this estimator is associated with, pointer makes
  // comparison fast
  const std::shared_ptr<const CSGSurface> surface;
};

/// @brief Collects multiple Estimator objects
class EstimatorSet {
  /// @brief EstimatorSetProxy objects will require access to many private
  ///        members of EstimatorSet
  friend EstimatorSetProxy;

public:
  /// @brief Constructs an EstimatorSet from an XML document
  /// @param estimators_node `estimators` node of an XML document
  /// @param world Used to reference CSGSurface or Cell objects for estimators
  /// @param perturbations Used to associate Estimator objects which have a
  ///        `sensitivities` node with a Perturbation
  /// @param total_weight Used for normalizing Estimator scores
  EstimatorSet(
      const pugi::xml_node& estimators_node, const World& world,
      const PerturbationSet& perturbations, const Real total_weight);
  /// @brief Copy constructor; deep copies each Estimator in other
  EstimatorSet(const EstimatorSet& other) noexcept;
  /// @brief Return an EstimatorSet::Proxy
  EstimatorSetProxy GetProxy() const noexcept;
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

/// @brief A lightweight class for accumulating scores produced by a single
///        history.
/// @details Constructing an EstimatorSet at the beginning of each history may
///          be expensive since it would require allocating a bunch of zero
///          scores
class EstimatorSetProxy {
public:
  /// @brief Construct a Proxy for the given EstimatorSet
  EstimatorSetProxy(const EstimatorSet& init) noexcept;
  /// @brief Score to each Estimator::Proxy
  void Score(const Particle& p) noexcept;
  /// @brief Calls Estimator::Proxy::CommitHistory() on each Estimator::Proxy
  void CommitHistory() noexcept;

private:
  // contains all estimator proxies
  std::vector<EstimatorProxy> estimator_proxies;
};
