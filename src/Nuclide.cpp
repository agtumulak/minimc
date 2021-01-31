#include "Nuclide.hpp"

// Nuclide

//// public

Nuclide::Nuclide(const pugi::xml_node& root, const std::string& nuclide_name)
    : name{nuclide_name}, xs{NuclearData::Create(root, nuclide_name)} {}

NuclearData::CrossSection Nuclide::GetTotal(const Particle& p) const noexcept {
  return xs.at(p.type)->GetTotal(p);
}

void Nuclide::Scatter(RNG& rng, Particle& p) const {
  xs.at(p.type)->Scatter(rng, p);
}

std::vector<Particle> Nuclide::Fission(RNG& rng, Particle& p) const noexcept {
  return xs.at(p.type)->Fission(rng, p);
}

NuclearData::Reaction
Nuclide::SampleReaction(RNG& rng, const Particle& p) const noexcept {
  return xs.at(p.type)->SampleReaction(rng, p);
}
