#include "Nuclide.hpp"

#include "History.hpp"
#include "State.hpp"
#include "pugixml.hpp"

#include <map>

// Nuclide

//// public

Nuclide::Nuclide(const pugi::xml_node& nuclide_node)
    : name{nuclide_node.attribute("name").as_string()}, xs{Interaction::Create(
                                                            nuclide_node)} {}

MicroscopicCrossSection Nuclide::GetMajorant(const State& s) const noexcept {
  return xs.at(s.particle)->GetMajorant(s);
}

MicroscopicCrossSection Nuclide::GetTotal(const State& s) const noexcept {
  return xs.at(s.particle)->GetTotal(s);
}

void Nuclide::Interact(State& s) const noexcept {
  xs.at(s.particle)->Interact(s);
}

void Nuclide::AppendSecondaries(
    std::vector<History> bank, const State& state) const noexcept {}
