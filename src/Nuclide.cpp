#include "Nuclide.hpp"

#include "InteractionDelegate.hpp"
#include "pugixml.hpp"

#include <cassert>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

// Nuclide

//// public

Nuclide::Nuclide(const pugi::xml_node& nuclide_node)
    : name{nuclide_node.attribute("name").as_string()},
      awr{nuclide_node.attribute("awr").as_double()}, xs{[this,
                                                          &nuclide_node]() {
        std::map<Particle::Type, std::unique_ptr<const InteractionDelegate>> xs;
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
                nuclide_node.path() + ": \"" + particle_name +
                "\" node not found");
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
                std::make_unique<const Continuous>(particle_node, *this));
          }
          else {
            assert(false); // this should have been caught by the validator
          }
        }
        return xs;
      }()} {}

MicroscopicCrossSection
Nuclide::GetCellMajorant(const Particle& p) const noexcept {
  return xs.at(p.type)->GetCellMajorant(p);
}

MicroscopicCrossSection Nuclide::GetTotal(const Particle& p) const noexcept {
  return xs.at(p.type)->GetTotal(p);
}

void Nuclide::Interact(
    Particle& p,
    std::vector<Estimator::Proxy>& estimator_proxies) const noexcept {
  xs.at(p.type)->Interact(p, estimator_proxies);
}

const std::optional<ThermalScattering>& Nuclide::GetTNSL() const {
  return xs.at(Particle::Type::neutron)->GetTNSL();
}
