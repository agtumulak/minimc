#include "Sensitivity.hpp"

#include "Estimator.hpp"
#include "Particle.hpp"
#include "pugixml.hpp"

#include <cassert>

// Sensitivity

//// public

std::unique_ptr<Sensitivity>
Sensitivity::Create(const pugi::xml_node& sensitivity_node) noexcept {
  switch (ToPerturbation(sensitivity_node.attribute("perturb").as_string())) {
  case Perturbation::total:
    break;
  case Perturbation::capture:
    break;
  case Perturbation::scatter:
    break;
  }

  if (perturbation_type == "sensitivity") {
    const std::string sensitivity_type =
        sensitivity_node.attribute("perturb").as_string();
    if (sensitivity_type == "total") {
      return std::make_unique<TotalCrossSectionPerturbation>();
    }
    else {
      assert(false); // this should have been caught by the validator
      return {};
    }
  }
  else {
    assert(false); // this should have been caught by the validator
    return {};
  }
}

Sensitivity::~Sensitivity() noexcept {}

// NoSensitivity

//// public

Real NoSensitivity::GetDirectEffect(
    const Particle&, const Estimator&) const noexcept {
  return 1;
}

// TotalCrossSectionPerturbation

//// public

Real TotalCrossSectionPerturbation::GetDirectEffect(
    const Particle& p, const Estimator& e) const noexcept {
  // get direct effect
  e.GetDirectEffect(p, *this);
  // get indirect effect
}
