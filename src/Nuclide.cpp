#include "Nuclide.hpp"

#include "Particle.hpp"
#include "Reaction.hpp"

// Nuclide

//// public

Nuclide::Nuclide(const pugi::xml_node& nuclide_node)
    : name{nuclide_node.attribute("name").as_string()}, xs{Interaction::Create(
                                                            nuclide_node)} {}

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

void Nuclide::Interact(Particle& p) const noexcept {
  xs.at(p.type)->Interact(p);
}
