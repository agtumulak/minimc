#include "Perturbation.hpp"

#include "Perturbation/IndirectEffect/IndirectEffect.hpp"
#include "Perturbation/Sensitivity/Sensitivity.hpp"
#include "World.hpp"
#include "pugixml.hpp"

#include <cassert>

using namespace Perturbation;

// Interface

//// public

std::unique_ptr<Interface> Interface::Create(
    const pugi::xml_node& perturbation_node, const World& world) noexcept {
  const std::string perturbation_type = perturbation_node.name();
  if (perturbation_type == "totalxs") {
    return std::make_unique<TotalCrossSection>(perturbation_node, world);
  }
  else {
    assert(false); // this should have been caught by the validator
    return {};
  }
}

Interface::~Interface() noexcept {}

//// protected

Interface::Interface(
    const pugi::xml_node& perturbation_node,
    const size_t n_perturbations) noexcept
    : name{perturbation_node.attribute("name").as_string()},
      n_perturbations{n_perturbations} {}

// TotalCrossSection

//// public

TotalCrossSection::TotalCrossSection(
    const pugi::xml_node& totalxs_node, const World& world)
    : Interface{totalxs_node, 1},
      nuclide{world.FindNuclideByName(
          totalxs_node.attribute("nuclide").as_string())} {}

std::unique_ptr<Sensitivity::Interface> TotalCrossSection::CreateSensitivity(
    const Estimator::Interface& estimator) const noexcept {
  return std::make_unique<Sensitivity::TotalCrossSection>(estimator, *this);
}

std::unique_ptr<IndirectEffect::Interface>
TotalCrossSection::CreateIndirectEffect() const noexcept {
  return std::make_unique<IndirectEffect::TotalCrossSection>(*this);
};
