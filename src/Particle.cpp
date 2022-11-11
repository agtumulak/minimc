#include "Particle.hpp"

#include "Constants.hpp"

#include <cassert>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

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
    return {};
  };
}

Particle::Particle(
    const std::unordered_map<const Perturbation*, Real>& indirect_effects,
    const Point& position, const Direction& direction, const Energy& energy,
    const Cell* cell, const Type type, RNG::result_type seed) noexcept
    : indirect_effects{indirect_effects}, position{position},
      direction{direction}, energy{energy}, cell{cell}, type{type}, rng{seed} {}

void Particle::Scatter(const Real& mu, const Energy& e) noexcept {
  // update direction and energy
  const Real phi = 2 * constants::pi * std::uniform_real_distribution{}(rng);
  direction = Direction{direction, mu, phi};
  energy = e;
}

const Point& Particle::GetPosition() const noexcept { return position; };

void Particle::SetPosition(const Point& p) noexcept { position = p; }

const Direction& Particle::GetDirection() const noexcept { return direction; };

void Particle::SetDirection(const Direction& d) noexcept { direction = d; }

void Particle::SetDirectionIsotropic() noexcept { direction = {rng}; }

const Energy& Particle::GetEnergy() const noexcept { return energy; };

void Particle::SetEnergy(const Energy& e) noexcept { energy = e; };

Particle::Type Particle::GetType() const noexcept { return type; }

const Cell& Particle::GetCell() const {
  if (cell) {
    return *cell;
  }
  throw std::runtime_error("Particle does not belong to a Cell");
};

void Particle::SetCell(const Cell& c) noexcept { cell = &c; };

void Particle::BankSecondaries(
    const Direction& direction, const Energy& energy) noexcept {
  secondaries.emplace_back(
      indirect_effects, position, direction, energy, cell, Type::neutron,
      rng());
}

void Particle::MoveSecondariesTo(Bank& bank) noexcept {
  bank += std::move(secondaries);
}

Real Particle::Sample() noexcept {
  return std::uniform_real_distribution{}(rng);
}

bool Particle::IsAlive() const noexcept {
  return event != Event::capture && event != Event::leak &&
         event != Event::fission;
}
