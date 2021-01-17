#include "Source.hpp"

#include "World.hpp"

#include <sstream>

// Source

//// public

Source::Source(const pugi::xml_node& root)
    : particle_type{CreateParticleType(root)}, initial_energy{
                                                   CreateDefaultEnergy(root)} {}

Particle Source::Sample(std::minstd_rand& rng, const World& w) const noexcept {
  Particle p{initial_energy, particle_type};
  p.SetDirectionIsotropic(rng);
  p.SetCell(w.FindCellContaining(p.GetPosition()));
  return p;
}

//// private

Particle::Type Source::CreateParticleType(const pugi::xml_node& root) noexcept {
  std::string particle_name;
  std::stringstream particle_name_list{
      root.child("general").child("particles").child_value()};
  particle_name_list >> particle_name;
  return Particle::ToType(particle_name);
}

Energy Source::CreateDefaultEnergy(const pugi::xml_node& root) noexcept {
  const std::string energy_type{root.child("nuclides").first_child().name()};
  if (energy_type == "multigroup") {
    return Group{1};
  }
  else if (energy_type == "continuous") {
    return ContinuousEnergy{1e-6};
  }
  else {
    assert(false); // this should hae been caught by the validator
  }
}
