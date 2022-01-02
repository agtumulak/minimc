#include "Sensitivity.hpp"

#include "Estimator.hpp"
#include "Particle.hpp"
#include "pugixml.hpp"

#include <string>

// Sensitivity

//// public

std::unique_ptr<Sensitivity> Sensitivity::Create(
    const pugi::xml_node& sensitivity_perturbation_node,
    const PerturbationSet& perturbations, const Estimator& estimator) {
  // find matching Perturbation
  const auto& matched_perturbation = perturbations.FindPerturbationByName(
      sensitivity_perturbation_node.attribute("name").as_string());
  return matched_perturbation.CreateSensitivity(
      sensitivity_perturbation_node, estimator);
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

// CurrentTotalCrossSectionSensitivity

//// public

CurrentTotalCrossSectionSensitivity::CurrentTotalCrossSectionSensitivity(
    const Estimator& estimator, const Perturbation& perturbation) noexcept
    : Sensitivity{estimator, perturbation} {}

std::unique_ptr<Sensitivity>
CurrentTotalCrossSectionSensitivity::Clone() const noexcept {
  return std::make_unique<CurrentTotalCrossSectionSensitivity>(*this);
}

Real CurrentTotalCrossSectionSensitivity::GetScore(
    const Particle& p) const noexcept {
  return p.GetIndirectEffect(&perturbation) * estimator.GetScore(p);
}
