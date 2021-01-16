#include "Nuclide.hpp"

// Nuclide

//// public

Nuclide::Nuclide(const pugi::xml_node& root, const std::string& nuclide_name)
    : name{nuclide_name}, xs{NuclearData::Create(root, nuclide_name)} {}

NuclearData::CrossSection Nuclide::GetTotal(const Particle& p) const noexcept {
  return xs.at(p.type)->GetTotal(p);
}

void Nuclide::Scatter(std::minstd_rand& rng, Particle& p) const {
  xs.at(p.type)->Scatter(rng, p);
}

NuclearData::Reaction Nuclide::SampleReaction(
    std::minstd_rand& rng, const Particle& p) const noexcept {
  return xs.at(p.type)->SampleReaction(rng, p);
}
