#include "Sensitivity.hpp"

#include "Estimator.hpp"
#include "Perturbation.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>

// Sensitivity

//// public

std::unique_ptr<Sensitivity> Sensitivity::Create(
    const pugi::xml_node& sensitivity_perturbation_node,
    const std::vector<std::unique_ptr<const Perturbation>>& perturbations,
    const Estimator& estimator) {
  // find matching Perturbation
  const std::string& name =
      sensitivity_perturbation_node.attribute("name").as_string();
  const auto perturbation_it = std::find_if(
      perturbations.cbegin(), perturbations.cend(),
      [&name](const auto& perturbation_ptr) {
        return perturbation_ptr->name == name;
      });
  if (perturbation_it == perturbations.cend()) {
    throw std::runtime_error(
        "Perturbation \"" + name + "\" not found. Must be one of: [" +
        std::accumulate(
            perturbations.cbegin(), perturbations.cend(), std::string{},
            [](const auto& accumulated, const auto& perturbation_ptr) noexcept {
              return accumulated + "\"" + perturbation_ptr->name + "\", ";
            }) +
        "]");
  }
  return (*perturbation_it)
      ->CreateSensitivity(sensitivity_perturbation_node, estimator);
}

Sensitivity::~Sensitivity() noexcept {}

std::string Sensitivity::to_string(const Real total_weight) const noexcept {
  return name + "\n" + std::string(name.size(), '=') + "\n\n" +
         GetScoreAsString(total_weight) + "\n";
}

Sensitivity& Sensitivity::operator+=(const Sensitivity& other) noexcept {
  Scorable::operator+=(other);
  return *this;
}

//// protected

Sensitivity::Sensitivity(
    const Estimator& estimator, const Perturbation& perturbation) noexcept
    : Scorable{estimator, perturbation}, estimator{estimator},
      perturbation{perturbation} {}

// SensitivityProxy

//// public

SensitivityProxy::SensitivityProxy(Sensitivity& original) noexcept
    : original{original} {};

void SensitivityProxy::CommitHistory() const noexcept {
  for (const auto& [index, score] : pending_scores) {
    original.AddScore(index, score);
  }
}

// CurrentTotalCrossSectionSensitivity

//// public

CurrentTotalCrossSectionSensitivity::CurrentTotalCrossSectionSensitivity(
    const Estimator& estimator, const Perturbation& perturbation) noexcept
    : Sensitivity{estimator, perturbation} {}

std::unique_ptr<Sensitivity>
CurrentTotalCrossSectionSensitivity::Clone() const noexcept {
  return std::make_unique<CurrentTotalCrossSectionSensitivity>(*this);
}
