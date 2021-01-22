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

Particle::Particle(
    const Point& position, const Direction& direction, const Energy& energy,
    const Type type) noexcept
    : position{position}, direction{direction}, energy{energy}, type{type} {}

bool Particle::IsAlive() const noexcept { return alive; }

void Particle::Kill() noexcept { alive = false; }

void Particle::Stream(const Real distance) noexcept {
  position += direction * distance;
}

const Point& Particle::GetPosition() const noexcept { return position; };

const Direction& Particle::GetDirection() const noexcept { return direction; };

void Particle::SetDirectionIsotropic(RNG& rng) noexcept {
  direction = Direction::CreateIsotropic(rng);
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
