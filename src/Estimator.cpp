#include "Estimator.hpp"

#include "Bins.hpp"
#include "Particle.hpp"
#include "World.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <ostream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>

// Estimator

std::ostream& operator<<(std::ostream& os, const Estimator& e) noexcept {
  os << "cosine: " << *e.cosine << std::endl;
  os << "energy: " << *e.energy << std::endl;
  for (const auto& score : e.scores) {
    os << score << " ";
  }
  return os;
}

//// public

std::unique_ptr<Estimator> Estimator::Create(
    const pugi::xml_node& estimator_node, const World& world) noexcept {
  const std::string estimator_type = estimator_node.name();
  if (estimator_type == "current") {
    return std::make_unique<CurrentEstimator>(estimator_node, world);
  }
  else {
    assert(false); // this should have been caught by the validator
  }
}

Estimator::~Estimator() noexcept {}

const std::vector<Real>& Estimator::GetScores() const noexcept {
  return scores;
}

Estimator& Estimator::Normalize(Real total_weight) noexcept {
  std::transform(
      scores.cbegin(), scores.cend(), scores.begin(),
      std::bind(std::divides<Real>(), std::placeholders::_1, total_weight));
  return *this;
}

Estimator& Estimator::operator+=(const Estimator& other) noexcept {
  std::transform(
      scores.cbegin(), scores.cend(), other.scores.cbegin(), scores.begin(),
      std::plus<Real>());
  return *this;
}

//// protected

Estimator::Estimator(const pugi::xml_node& estimator_node) noexcept
    : direction{CreateDirection(estimator_node.child("cosine"))},
      cosine{Bins::Create(estimator_node.child("cosine").first_child())},
      energy{Bins::Create(estimator_node.child("energy").first_child())},
      strides{ComputeStrides(*cosine, *energy)},
      scores(cosine->size() * strides.front(), 0) {}

Real& Estimator::GetScore(const Particle& p) noexcept {
  // get cosine bin
  const auto c_i =
      direction ? cosine->GetIndex(direction.value().Dot(p.GetDirection())) : 0;
  // get energy bin
  const auto e_i = energy->GetIndex(std::visit(VisitEnergy(), p.GetEnergy()));
  // get flattened index
  return scores[strides[0] * c_i + e_i];
}

//// private

std::optional<Direction>
Estimator::CreateDirection(const pugi::xml_node& cosine_node) noexcept {
  return cosine_node
             ? std::optional<Direction>(
                   std::in_place,
                   cosine_node.attribute("u").as_double(),
                   cosine_node.attribute("v").as_double(),
                   cosine_node.attribute("w").as_double())
             : std::nullopt;
}

// CurrentEstimator

//// public

CurrentEstimator::CurrentEstimator(
    const pugi::xml_node& current_estimator_node, const World& world)
    : Estimator{current_estimator_node},
      surface{world.FindSurfaceByName(
          current_estimator_node.attribute("surface").as_string())} {}

std::unique_ptr<Estimator> CurrentEstimator::Clone() const noexcept {
  return std::make_unique<CurrentEstimator>(*this);
}

void CurrentEstimator::Score(const Particle& p) noexcept {
  if (p.current_surface == surface &&
      (p.event == Particle::Event::surface_cross ||
       p.event == Particle::Event::leak)) {
    GetScore(p) += 1;
  }
}

// EstimatorSet

std::ostream& operator<<(std::ostream& os, const EstimatorSet& e) noexcept {
  for (const auto& [name, estimator] : e.estimators) {
    os << std::endl
       << name << std::endl
       << std::string(name.length(), '-') << std::endl
       << *estimator;
  }
  return os;
}

//// public


EstimatorSet::EstimatorSet(const EstimatorSet& other) noexcept
    : estimators{ConstructAllEstimators(other)} {}

EstimatorSet::EstimatorSet(
    const pugi::xml_node& estimators_node, const World& world)
    : estimators{ConstructAllEstimators(estimators_node, world)} {}

const Estimator& EstimatorSet::GetEstimator(const std::string& name) const {
  return *estimators.at(name);
}

void EstimatorSet::Score(const Particle& p) noexcept {
  std::for_each(
      estimators.begin(), estimators.end(), [&p](auto& name_estimator_pair) {
        auto& [name, estimator] = name_estimator_pair;
        estimator->Score(p);
      });
}

EstimatorSet& EstimatorSet::Normalize(Real total_weight) noexcept {
  std::for_each(
      estimators.begin(), estimators.end(),
      [total_weight](auto& name_estimator_pair) {
        auto& [name, estimator] = name_estimator_pair;
        estimator->Normalize(total_weight);
      });
  return *this;
}

EstimatorSet& EstimatorSet::operator+=(const EstimatorSet& other) noexcept {
  std::for_each(
      estimators.begin(), estimators.end(),
      [&other](auto& name_estimator_pair) {
        auto& [name, estimator] = name_estimator_pair;
        *estimator += *(other.estimators.at(name));
      });
  return *this;
}

//// private

std::map<std::string, std::unique_ptr<Estimator>>
EstimatorSet::ConstructAllEstimators(const EstimatorSet& other) noexcept {
  std::map<std::string, std::unique_ptr<Estimator>> result;
  for (const auto& [name, estimator] : other.estimators) {
    result[name] = estimator->Clone();
  }
  return result;
}

std::map<std::string, std::unique_ptr<Estimator>>
EstimatorSet::ConstructAllEstimators(
    const pugi::xml_node& estimators_node, const World& world) {
  std::map<std::string, std::unique_ptr<Estimator>> result;
  for (const auto& estimator_node : estimators_node) {
    // check if an Estimator with this name was already constructed
    if (const std::string estimator_name =
            estimator_node.attribute("name").as_string();
        result.find(estimator_name) != result.cend()) {
      throw std::runtime_error("Duplicate estimator name: " + estimator_name);
    }
    else {
      result[estimator_name] = Estimator::Create(estimator_node, world);
    }
  }
  return result;
}
