#include "Nuclide.hpp"

#include "Particle.hpp"
#include "pugixml.hpp"

#include <map>

// Nuclide

//// public

Nuclide::Nuclide(const pugi::xml_node& nuclide_node)
    : name{nuclide_node.attribute("name").as_string()}, xs{Interaction::Create(
                                                            nuclide_node)} {}

MicroscopicCrossSection Nuclide::GetMajorant(const Particle& p) const noexcept {
  return xs.at(p.type)->GetMajorant(p);
}

MicroscopicCrossSection Nuclide::GetTotal(const Particle& p) const noexcept {
  return xs.at(p.type)->GetTotal(p);
}

void Nuclide::Interact(Particle& p) const noexcept {
  xs.at(p.type)->Interact(p);
}
