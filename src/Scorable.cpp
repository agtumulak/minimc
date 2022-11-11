#include "Scorable.hpp"

#include "Bins.hpp"
#include "Perturbation.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <sstream>
#include <string>

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

void Scorable::AddScore(const size_t index, const Real score) noexcept {
  scores[index] += score;
  square_scores[index] += score * score;
}

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
  std::stringstream sstream;
  sstream << "mean\n----\n";
  sstream << std::scientific;
  for (const auto& score : scores) { // wow so clean
    sstream << score / total_weight << ", ";
  }
  sstream << "\n\n";
  // gross
  sstream << "std dev\n-------\n";
  for (size_t i = 0; i < scores.size(); i++) { // gross
    const auto& score = scores[i];
    const auto& square_score = square_scores[i];
    sstream << std::sqrt(square_score - score * score / total_weight) /
                   total_weight
            << ", ";
  }
  sstream << "\n\n";
  return sstream.str();
}
