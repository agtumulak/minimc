#include "Particle.hpp"

#include <stdexcept>

// Particle

//// public

Particle::Type Particle::ToType(const std::string& name) noexcept {
  if (name == "neutron") {
    return Particle::Type::neutron;
  }
  else if (name == "photon") {
    return Particle::Type::photon;
  }
  else {
    assert(false); // this should have been caught by the validator
  };
}

Particle::Particle() noexcept {}

Particle::Particle(const Energy& energy, const Type type) noexcept
    : energy{energy}, type{type} {};

Particle::operator bool() const noexcept { return alive; }

const Point& Particle::GetPosition() const noexcept { return position; };

const Point& Particle::GetDirection() const noexcept { return direction; };

void Particle::SetDirectionIsotropic(std::minstd_rand& rng) noexcept {
  direction.SetIsotropic(rng);
}

const Energy& Particle::GetEnergy() const noexcept { return energy; };

void Particle::SetEnergy(const Energy& e) noexcept { energy = e; };

const Cell& Particle::GetCell() const {
  if (cell) {
    return *cell;
  }
  throw std::runtime_error("Particle does not belong to a Cell");
}

void Particle::SetCell(const Cell& c) noexcept { cell = &c; };
