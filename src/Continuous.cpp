#include "Continuous.hpp"

#include <fstream>
#include <limits>
#include <string>
#include <variant>

// Continuous

//// public

Continuous::Continuous(const pugi::xml_node& particle_node)
    : reactions{CreateReactions(particle_node)},
      total{particle_node.child("total").attribute("file").as_string()} {}

Continuous::CrossSection
Continuous::GetTotal(const Particle& p) const noexcept {
  return total.at(std::get<ContinuousEnergy>(p.GetEnergy()));
}

void Continuous::Scatter(std::minstd_rand& rng, Particle& p) const noexcept {
  p.SetDirectionIsotropic(rng);
  return;
}

NuclearData::Reaction Continuous::SampleReaction(
    std::minstd_rand& rng, const Particle& p) const noexcept {
  const CrossSection threshold =
      std::uniform_real_distribution{}(rng)*GetTotal(p);
  CrossSection accumulated{0};
  for (const auto& [reaction, xs] : reactions) {
    accumulated += xs.at(std::get<ContinuousEnergy>(p.GetEnergy()));
    if (accumulated > threshold) {
      return reaction;
    }
  }
  // If no reaction found, resample tail-recursively
  return SampleReaction(rng, p);
}

// Continuous::OneDimensional

//// public

Continuous::OneDimensional::OneDimensional(
    const std::filesystem::path& datapath) {
  std::ifstream datafile{datapath};
  if (!datafile) {
    throw std::runtime_error("File not found: " + datapath.string());
  }
  for (size_t i = 1; i <= 3; i++) {
    // first three lines are header
    datafile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  ContinuousEnergy energy;
  NuclearData::CrossSection xs;
  while (datafile >> energy) {
    datafile.ignore(
        std::numeric_limits<std::streamsize>::max(), ';'); // delimiter
    datafile >> xs;
    elements[energy] = xs;
  }
}

const Real&
Continuous::OneDimensional::at(const ContinuousEnergy e) const noexcept {
  // TODO: Interpolate and handle edge cases
  return (*elements.upper_bound(e)).second;
}

//  Continuous

//// private

Continuous::ReactionsMap
Continuous::CreateReactions(const pugi::xml_node& particle_node) {
  Continuous::ReactionsMap reactions;
  for (const auto& reaction_node : particle_node) {
    const std::string reaction_name = reaction_node.name();
    if (reaction_name == "total") {
      continue; // skip total cross section
    }
    reactions.emplace(
        NuclearData::ToReaction(reaction_name),
        OneDimensional{reaction_node.attribute("file").as_string()});
  }
  return reactions;
}
