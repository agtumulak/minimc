#include "Estimator.hpp"

#include "Bins.hpp"
#include "Particle.hpp"
#include "World.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <type_traits>
#include <utility>
#include <variant>

// Estimator

//// public

std::unique_ptr<Estimator> Estimator::Create(
    const pugi::xml_node& estimator_node, const World& world) noexcept {
  const std::string estimator_type = estimator_node.name();
  if (estimator_type == "current") {
    return std::make_unique<CurrentEstimator>(estimator_node, world);
  }
  else {
    assert(false); // this should have been caught by the validator
    return {};
  }
}

Estimator::~Estimator() noexcept {}

std::ostringstream
Estimator::GetPrintable(const Real total_weight) const noexcept {
  std::ostringstream ostream;
  ostream << "cosine\n------\n" << *cosine << "\n" << std::endl;
  ostream << "energy\n------\n" << *energy << "\n" << std::endl;
  ostream << "scores\n------\n";
  // normalize each score by total_weight
  std::vector<Real> estimates(scores.size(), 0);
  std::transform(
      scores.cbegin(), scores.cend(), estimates.begin(),
      std::bind(std::divides<Real>(), std::placeholders::_1, total_weight));
  for (const auto& estimate : estimates) {
    ostream << estimate << ", ";
  }
  // compute standard deviation of the mean
  std::vector<Real> stddevs(scores.size(), 0);
  std::transform(
      scores.cbegin(), scores.cend(), square_scores.cbegin(), stddevs.begin(),
      [total_weight](const Real& score, const Real& square_score) {
        return std::sqrt(square_score - score * score / total_weight) /
               total_weight;
      });
  ostream << "\n" << std::endl;
  ostream << "standard deviations\n-------------------\n";
  for (const auto& stddev : stddevs) {
    ostream << stddev << ", ";
  }
  ostream << "\n" << std::endl;
  return ostream;
}

Estimator& Estimator::operator+=(const Estimator& other) noexcept {
  // add scores
  std::transform(
      scores.cbegin(), scores.cend(), other.scores.cbegin(), scores.begin(),
      std::plus<Real>());
  // add square scores
  std::transform(
      square_scores.cbegin(), square_scores.cend(),
      other.square_scores.cbegin(), square_scores.begin(), std::plus<Real>());
  return *this;
}

//// protected

Estimator::Estimator(const pugi::xml_node& estimator_node) noexcept
    : name{estimator_node.attribute("name").as_string()},
      direction{CreateDirection(estimator_node.child("cosine"))},
      cosine{Bins::Create(estimator_node.child("cosine").first_child())},
      energy{Bins::Create(estimator_node.child("energy").first_child())},
      strides{ComputeStrides(*cosine, *energy)},
      scores(cosine->size() * strides.front(), 0),
      square_scores(cosine->size() * strides.front(), 0) {}

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

size_t Estimator::GetIndex(const Particle& p) const noexcept {
  // get cosine bin
  const auto c_i =
      direction ? cosine->GetIndex(direction.value().Dot(p.GetDirection())) : 0;
  // get energy bin
  const auto e_i = energy->GetIndex(std::visit(VisitEnergy(), p.GetEnergy()));
  // get flattened index
  return strides[0] * c_i + e_i;
}

// Estimator::Proxy

//// public

Estimator::Proxy::Proxy(Estimator& init) noexcept : original{init} {}

void Estimator::Proxy::Score(const Particle& p) noexcept {
  // determine score this Particle would produce
  const auto& score = original.GetScore(p);
  // skip scoring if score would have been zero
  if (score == 0) {
    return;
  }
  // determine which bin this Particle would score to
  const auto& index = original.GetIndex(p);
  // check if the Particle's bin has already been scored to
  if (const auto it = pending_scores.find(index); it != pending_scores.cend()) {
    // add score to existing index
    it->second += score;
  }
  else {
    // insert index and initialize score
    pending_scores.insert(std::make_pair(index, score));
  }
}

void Estimator::Proxy::CommitHistory() const noexcept {
  for (const auto& [index, score] : pending_scores) {
    original.scores[index] += score;
    original.square_scores[index] += score * score;
  }
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

Real CurrentEstimator::GetScore(const Particle& p) noexcept {
  if (p.current_surface == surface &&
      (p.event == Particle::Event::surface_cross ||
       p.event == Particle::Event::leak)) {
    return 1;
  }
  else {
    return 0;
  }
}

// EstimatorSet

std::ostream& operator<<(std::ostream& os, const EstimatorSet& e) noexcept {
  for (const auto& estimator : e.estimators) {
    os << std::endl
       << estimator->name << std::endl
       << std::string(estimator->name.length(), '-') << std::endl
       << estimator->GetPrintable(e.total_weight).str() << std::endl;
  }
  return os;
}

//// public

EstimatorSet::EstimatorSet(
    const pugi::xml_node& estimators_node, const World& world,
    const Real total_weight)
    : // IIFE
      estimators{[&estimators_node, &world]() {
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

EstimatorSet::Proxy EstimatorSet::GetProxy() const noexcept {
  return Proxy(*this);
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

// EstimatorSet::Proxy

//// public

EstimatorSet::Proxy::Proxy(const EstimatorSet& init) noexcept
    : // IIFE
      estimator_proxies{[&init]() {
        std::vector<Estimator::Proxy> result;
        for (const auto& estimator_ptr : init.estimators) {
          result.emplace_back(*estimator_ptr);
        }
        return result;
      }()} {}

void EstimatorSet::Proxy::Score(const Particle& p) noexcept {
  std::for_each(
      estimator_proxies.begin(), estimator_proxies.end(),
      [&p](auto& estimator_proxy) { estimator_proxy.Score(p); });
}

void EstimatorSet::Proxy::CommitHistory() const noexcept {
  std::for_each(
      estimator_proxies.begin(), estimator_proxies.end(),
      [](auto& proxy) { proxy.CommitHistory(); });
}

