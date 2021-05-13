#include "NuclearData.hpp"

#include "Continuous.hpp"
#include "Multigroup.hpp"

#include <numeric>
#include <sstream>
#include <stdexcept>

// NuclearData

//// public

NuclearData::Map NuclearData::Create(
    const pugi::xml_node& root, const std::string& nuclide_name) {
  Map xs;
  const auto& energy_type_node{root.child("nuclides").first_child()};
  const auto& nuclide_node{
      energy_type_node.find_child_by_attribute("name", nuclide_name.c_str())};
  if (!nuclide_node) {
    throw std::runtime_error(
        "Nuclide node \"" + nuclide_name + "\" not found. Must be one of: [" +
        std::accumulate(
            energy_type_node.begin(), energy_type_node.end(), std::string{},
            [](const auto& accumulated, const auto& nuclide_node) noexcept {
              return accumulated + "\"" +
                     nuclide_node.attribute("name").as_string() + "\", ";
            }) +
        "]");
  }

  std::stringstream particle_name_list{
      root.child("general").child("particles").child_value()};
  std::string particle_name;
  while (particle_name_list >> particle_name) {
    const auto& particle_node{nuclide_node.child(particle_name.c_str())};
    if (!particle_node) {
      throw std::runtime_error(
          nuclide_node.path() + ": \"" + particle_name + "\" node not found");
    }
    const std::string energy_type{energy_type_node.name()};
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

NuclearData::~NuclearData() noexcept {}

