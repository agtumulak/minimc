#include "Scorable.hpp"

#include "Bins.hpp"
#include "Perturbation.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>

// Scorable

//// public

Scorable::Scorable(
    const std::string& name, const pugi::xml_node& bins_node) noexcept
    : name{name}, bins{std::make_shared<ParticleBins>(bins_node)},
      scores(bins->size(), 0), square_scores(bins->size(), 0) {}

Scorable::Scorable(
    const Scorable& estimator, const Perturbation& perturbation) noexcept
    : // The XML schema enforces that each Estimator has a unique name and each
      // Perturbation has a unique name. Concatenating the two will produce
      // a unique sensitivity name.
      name{estimator.name + "::" + perturbation.name}, bins{estimator.bins},
      scores(bins->size(), 0), square_scores(bins->size(), 0) {}

Scorable::~Scorable() noexcept {}

Real Scorable::GetScore(
    const size_t index, const Real total_weight) const noexcept {
  return scores.at(index) / total_weight;
}

Scorable& Scorable::operator+=(const Scorable& other) noexcept {
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

std::string Scorable::GetScoreAsString(const Real total_weight) const noexcept {
  std::string result;
  result += "mean\n----\n";
  for (const auto& score : scores) { // wow so clean
    result += std::to_string(score / total_weight) + ", ";
  }
  result += "\n\n";
  // gross
  result += "std dev\n-------\n";
  for (size_t i = 0; i < scores.size(); i++) { // gross
    const auto& score = scores[i];
    const auto& square_score = square_scores[i];
    result += std::to_string(
                  std::sqrt(square_score - score * score / total_weight) /
                  total_weight) +
              ", ";
  }
  result += "\n\n";
  return result;
}

// ScorableProxy

//// public

ScorableProxy::ScorableProxy(Scorable& original) noexcept
    : original{original} {}

ScorableProxy::~ScorableProxy() noexcept {}

void ScorableProxy::Score(const Particle& p) noexcept {
  // determine score this Particle would produce
  const auto& score = original.GetScore(p);
  // skip scoring if score would have been zero
  if (score == 0) {
    return;
  }
  // determine which bin this Particle would score to
  const auto& index = original.bins->GetIndex(p);
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

void ScorableProxy::CommitHistory() const noexcept {
  for (const auto& [index, score] : pending_scores) {
    original.scores[index] += score;
    original.square_scores[index] += score * score;
  }
}
