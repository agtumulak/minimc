#include "Interaction.hpp"

#include "Continuous.hpp"
#include "Multigroup.hpp"

#include <cassert>
#include <sstream>
#include <stdexcept>
#include <string>

// Interaction

//// public

Interaction::Map Interaction::Create(const pugi::xml_node& nuclide_node) {
  Map xs;
  std::stringstream particle_name_list{
      nuclide_node.root()
          .select_node("minimc/general/particles")
          .node()
          .child_value()};
  std::string particle_name;
  while (particle_name_list >> particle_name) {
    const auto& particle_node{nuclide_node.child(particle_name.c_str())};
    if (!particle_node) {
      throw std::runtime_error(
          nuclide_node.path() + ": \"" + particle_name + "\" node not found");
    }
    const std::string energy_type{nuclide_node.parent().name()};
    if (energy_type == "multigroup") {
      xs.emplace(
          Particle::ToType(particle_name),
          std::make_unique<const Multigroup>(particle_node));
    }
    else if (energy_type == "continuous") {
      xs.emplace(
          Particle::ToType(particle_name),
          std::make_unique<const Continuous>(particle_node));
    }
    else {
      assert(false); // this should have been caught by the validator
    }
  }
  return xs;
}

//// protected

Interaction::~Interaction() noexcept {}
