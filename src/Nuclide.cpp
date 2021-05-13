#include "Nuclide.hpp"

#include "Particle.hpp"

// Nuclide

//// public

Nuclide::Nuclide(const pugi::xml_node& root, const std::string& nuclide_name)
    : name{nuclide_name}, xs{NuclearData::Create(root, nuclide_name)} {}

MicroscopicCrossSection Nuclide::GetTotal(const Particle& p) const noexcept {
  return xs.at(p.type)->GetTotal(p);
}

MicroscopicCrossSection
Nuclide::GetReaction(const Particle& p, const Reaction r) const noexcept {
  return xs.at(p.type)->GetReaction(p, r);
}

Real Nuclide::GetNuBar(const Particle& p) const noexcept {
  return xs.at(p.type)->GetNuBar(p);
}

void Nuclide::Scatter(RNG& rng, Particle& p) const {
  xs.at(p.type)->Scatter(rng, p);
}

std::vector<Particle> Nuclide::Fission(RNG& rng, Particle& p) const noexcept {
  return xs.at(p.type)->Fission(rng, p);
}
