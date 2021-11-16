#include "Particle.hpp"

#include <stdexcept>
#include <string>

Particle ToParticle(const std::string& name) {
  if (name == "neutron") {
    return Particle::neutron;
  }
  else if (name == "photon") {
    return Particle::photon;
  }
  else {
    throw std::runtime_error("Unrecognized particle name: " + name);
  };
}
