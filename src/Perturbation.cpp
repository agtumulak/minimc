#include "Perturbation.hpp"

#include "IndirectEffect.hpp"
#include "Sensitivity.hpp"
#include "World.hpp"
#include "pugixml.hpp"

#include <cassert>
#include <stdexcept>

// Perturbation

//// public

std::unique_ptr<const Perturbation> Perturbation::Create(
    const pugi::xml_node& perturbation_node, const World& world) noexcept {
  const std::string perturbation_type = perturbation_node.name();
  if (perturbation_type == "total") {
    return std::make_unique<TotalCrossSectionPerturbation>(
        perturbation_node, world);
  }
  else {
    assert(false); // this should have been caught by the validator
    return {};
  }
}

Perturbation::~Perturbation() noexcept {}

//// protected

Perturbation::Perturbation(const pugi::xml_node& perturbation_node) noexcept
    : name{perturbation_node.attribute("name").as_string()} {}

// TotalCrossSectionPerturbation

//// public

TotalCrossSectionPerturbation::TotalCrossSectionPerturbation(
    const pugi::xml_node& total_node, const World& world)
    : Perturbation{total_node},
      nuclide{world.FindNuclideByName(
          total_node.attribute("nuclide").as_string())} {}

std::unique_ptr<IndirectEffect>
TotalCrossSectionPerturbation::CreateIndirectEffect() const noexcept {
  return std::make_unique<TotalCrossSectionPerturbationIndirectEffect>(*this);
};

std::unique_ptr<Sensitivity> TotalCrossSectionPerturbation::CreateSensitivity(
    const pugi::xml_node& perturbation_node, const Estimator& estimator) const {
  const pugi::xml_node& estimator_node = perturbation_node.parent().parent();
  const std::string estimator_type = estimator_node.name();
  if (estimator_type == "current") {
    return std::make_unique<CurrentTotalCrossSectionSensitivity>(
        estimator, *this);
  }
  else {
    throw std::runtime_error(
        std::string{"Sensitivity to total cross section perturbation not "
                    "supported for estimator: "} +
        estimator_node.path());
  }
}
