#include "Estimator.hpp"

#include "Bins.hpp"
#include "World.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

// Estimator::Visitor

//// public

Estimator::Visitor::~Visitor() noexcept {}

// Estimator

//// public

std::unique_ptr<Estimator> Estimator::Create(
    const pugi::xml_node& estimator_node, const World& world,
    const std::vector<std::unique_ptr<const Perturbation>>& perturbations) {
  const std::string estimator_type = estimator_node.name();
  if (estimator_type == "current") {
    return std::make_unique<CurrentEstimator>(
        estimator_node, world, perturbations);
  }
  else {
    assert(false); // this should have been caught by the validator
    return {};
  }
}

Estimator::Estimator(const Estimator& other) noexcept
    : Scorable(other),
      // IIFE
      sensitivities{[&other]() noexcept {
        std::vector<std::unique_ptr<Sensitivity>> result;
        for (const auto& sensitivity : other.sensitivities) {
          result.push_back(sensitivity->Clone());
        }
        return result;
      }()} {}

Estimator::~Estimator() noexcept {}

std::string Estimator::to_string(const Real total_weight) const noexcept {
  std::string result;
  result += name + "\n" + std::string(name.size(), '=') + "\n\n";
  result += bins->to_string();
  result += GetScoreAsString(total_weight);
  for (const auto& sensitivity : sensitivities) {
    result += sensitivity->to_string(total_weight);
  }
  return result;
}

Estimator& Estimator::operator+=(const Estimator& other) noexcept {
  // add scores
  Scorable::operator+=(other);
  // add sensitivity scores
  for (auto& sensitivity : sensitivities) {
    const auto& matched_sensitivity = *std::find_if(
        other.sensitivities.cbegin(), other.sensitivities.cend(),
        [&sensitivity](const auto& other_it) {
          return sensitivity->name == other_it->name;
        });
    *sensitivity += *matched_sensitivity;
  }
  return *this;
}

//// protected

Estimator::Estimator(
    const pugi::xml_node& estimator_node,
    const std::vector<std::unique_ptr<const Perturbation>>&
        perturbations) noexcept
    : Scorable{estimator_node.attribute("name").as_string(), estimator_node.child("bins")},
      // IIFE
      sensitivities{[this, &estimator_node, &perturbations]() {
        std::vector<std::unique_ptr<Sensitivity>> result;
        for (const auto& perturbation_sensitivity_node :
             estimator_node.child("sensitivities")) {
          result.push_back(Sensitivity::Create(
              perturbation_sensitivity_node, perturbations, *this));
        }
        return result;
      }()} {}

// EstimatorProxy

//// public

EstimatorProxy::EstimatorProxy(Estimator& original) noexcept
    : // IIFE
      sensitivity_proxies{[&original]() noexcept {
        std::vector<SensitivityProxy> result;
        for (const auto& sensitivity : original.sensitivities) {
          result.emplace_back(*sensitivity);
        }
        return result;
      }()},
      original{original} {}

void EstimatorProxy::Visit(const Estimator::Visitor& visitor) noexcept {
  if (const auto result = original.Visit(visitor)) {
    const auto& [bin_index, score] = result.value();
    // check if the Particle's bin has already been scored to
    if (const auto it = pending_scores.find(bin_index);
        it != pending_scores.cend()) {
      // add score to existing index
      it->second += score;
    }
    else {
      // insert index and initialize score
      pending_scores.insert(std::make_pair(bin_index, score));
    }
  }
}

void EstimatorProxy::CommitHistory() const noexcept {
  for (const auto& [index, score] : pending_scores) {
    original.AddScore(index, score);
  }
  // commit child ScorableProxy objects for Sensitivity objects
  for (auto& sensitivity_proxy : sensitivity_proxies) {
    sensitivity_proxy.CommitHistory();
  }
}

// CurrentEstimator

//// public

CurrentEstimator::CurrentEstimator(
    const pugi::xml_node& current_estimator_node, const World& world,
    const std::vector<std::unique_ptr<const Perturbation>>& perturbations)
    : Estimator{current_estimator_node, perturbations},
      surface{world.FindSurfaceByName(
          current_estimator_node.attribute("surface").as_string())} {}

CurrentEstimator::CurrentEstimator(const CurrentEstimator& other) noexcept
    : Estimator(other), surface{other.surface} {}

Estimator::Visitor::T
CurrentEstimator::Visit(const Visitor& visitor) const noexcept {
  return visitor.Visit(*this);
}

std::unique_ptr<Estimator> CurrentEstimator::Clone() const noexcept {
  return std::make_unique<CurrentEstimator>(*this);
}

// EstimatorSet

//// public

EstimatorSet::EstimatorSet(
    const pugi::xml_node& estimators_node, const World& world,
    const std::vector<std::unique_ptr<const Perturbation>>& perturbations,
    const Real total_weight)
    : // IIFE
      estimators{[&estimators_node, &world, &perturbations]() {
        std::vector<std::unique_ptr<Estimator>> result;
        for (const auto& estimator_node : estimators_node) {
          result.push_back(
              Estimator::Create(estimator_node, world, perturbations));
        }
        return result;
      }()},
      total_weight{total_weight} {}

EstimatorSet::EstimatorSet(const EstimatorSet& other) noexcept
    : // IIFE
      estimators{[&other]() {
        std::vector<std::unique_ptr<Estimator>> result;
        for (const auto& estimator : other.estimators) {
          result.push_back(estimator->Clone());
        }
        return result;
      }()},
      total_weight{other.total_weight} {}

std::vector<EstimatorProxy> EstimatorSet::CreateProxies() const noexcept {
  std::vector<EstimatorProxy> result;
  for (const auto& estimator : estimators) {
    result.emplace_back(*estimator);
  }
  return result;
}

const Estimator&
EstimatorSet::FindEstimatorByName(const std::string& name) const {
  const auto estimator_it = std::find_if(
      estimators.cbegin(), estimators.cend(),
      [&name](const auto& estimator_it) { return estimator_it->name == name; });
  if (estimator_it == estimators.cend()) {
    throw std::runtime_error(
        "Estimator \"" + name + "\" notfound. Must be one of: [" +
        std::accumulate(
            estimators.cbegin(), estimators.cend(), std::string{},
            [](const auto& accumulated, const auto& estimator_ptr) noexcept {
              return accumulated + "\"" + estimator_ptr->name + "\", ";
            }) +
        "]");
  }
  else {
    return **estimator_it;
  }
}

std::string EstimatorSet::to_string() const noexcept {
  std::string result;
  for (const auto& estimator : estimators) {
    result += "\n" + estimator->to_string(total_weight) + "\n";
  }
  return result;
}

EstimatorSet& EstimatorSet::operator+=(const EstimatorSet& other) {
  std::for_each(
      estimators.begin(), estimators.end(), [&other](auto& estimator) {
        // find the Estimator in other with a matching name
        const auto matched_it = std::find_if(
            other.estimators.cbegin(), other.estimators.cend(),
            [&estimator](const auto& other_it) {
              return estimator->name == other_it->name;
            });
        // throw exception if no matching Estimator was found
        if (matched_it == other.estimators.cend()) {
          throw std::runtime_error("Estimator not found: " + estimator->name);
        }
        *estimator += **matched_it;
      });
  return *this;
}
