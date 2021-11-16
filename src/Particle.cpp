#include "Particle.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "Material.hpp"
#include "Nuclide.hpp"

#include <cassert>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>

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
    const Point& position, const Direction& direction, const Energy& energy,
    const Type type, RNG::result_type seed, const Cell* cell) noexcept
    : position{position}, direction{direction}, energy{energy}, cell{cell},
      rng{seed}, type{type} {}

void Particle::Stream(const Real distance) noexcept {
  position += direction * distance;
}

void Particle::Scatter(const Real& mu, const Energy& e) noexcept {
  const Real phi = 2 * constants::pi * std::uniform_real_distribution{}(rng);
  direction = Direction{direction, mu, phi};
  energy = e;
}

const Point& Particle::GetPosition() const noexcept { return position; };

const Direction& Particle::GetDirection() const noexcept { return direction; };

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
      position, direction, energy, Type::neutron, rng(), cell);
}

void Particle::MoveSecondariesToFrontOf(std::list<Particle>& bank) noexcept {
  bank.splice(bank.begin(), secondaries);
}

Real Particle::Sample() noexcept {
  return std::uniform_real_distribution{}(rng);
}

const Nuclide& Particle::SampleNuclide() noexcept {
  const MicroscopicCrossSection threshold =
      cell->material->GetMicroscopicTotal(*this) *
      std::uniform_real_distribution{}(rng);
  MicroscopicCrossSection accumulated = 0;
  for (const auto& [nuclide_ptr, afrac] : cell->material->afracs) {
    accumulated += afrac * nuclide_ptr->GetTotal(*this);
    if (accumulated > threshold) {
      return *nuclide_ptr;
    }
  }
  // This should never be reached since Material total cross section is
  // computed from constituent Nuclide total cross sections
  assert(false);
}

bool Particle::IsAlive() const noexcept {
  return event != Event::capture && event != Event::leak &&
         event != Event::fission;
}
