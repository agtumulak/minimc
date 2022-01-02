#include "Perturbation.hpp"

#include "Cell.hpp"
#include "Material.hpp"
#include "Particle.hpp"
#include "Sensitivity.hpp"
#include "TransportMethod.hpp"
#include "World.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

// Perturbation

//// public

std::unique_ptr<Perturbation> Perturbation::Create(
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

Real TotalCrossSectionPerturbation::Stream(
    const Particle& p, const Real distance) const noexcept {
  // if Particle is not streaming in a Material which contains the perturbed
  // Nuclide, there is no indirect effect
  const auto& atom_fractions = p.GetCell().material->afracs;
  if (atom_fractions.find(nuclide) == atom_fractions.cend()) {
    return 0;
  }
  return 1 / Particle::transport_method->GetCollisionProbabilityDensity(p) -
         distance;
}

Real TotalCrossSectionPerturbation::Scatter(
    const Particle&, const Real&, const Energy&) const noexcept {
  return 0;
}

// PerturbationSet

//// public

PerturbationSet::PerturbationSet(
    const pugi::xml_node& perturbations_node, const World& world) noexcept
    : perturbations{[&perturbations_node, &world]() {
        std::vector<std::unique_ptr<const Perturbation>> result;
        for (const auto& perturbation_node : perturbations_node) {
          result.push_back(Perturbation::Create(perturbation_node, world));
        }
        return result;
      }()} {}

const Perturbation&
PerturbationSet::FindPerturbationByName(const std::string& name) const {
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
  else {
    return **perturbation_it;
  }
}

std::map<const Perturbation*, Real>
PerturbationSet::GetIndirectEffects() const noexcept {
  std::map<const Perturbation*, Real> result;
  for (const auto& perturbation : perturbations) {
    result.insert(std::make_pair(perturbation.get(), 0));
  }
  return result;
}
