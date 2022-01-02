#include "Particle.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "Material.hpp"
#include "Nuclide.hpp"
#include "Perturbation.hpp"
#include "TransportMethod.hpp"

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

std::unique_ptr<const TransportMethod> Particle::transport_method;

Particle::Particle(
    const Point& position, const Direction& direction, const Energy& energy,
    const Type type, RNG::result_type seed, const Cell* cell) noexcept
    : position{position}, direction{direction}, energy{energy}, cell{cell},
      rng{seed}, type{type} {}

Bank Particle::Transport(EstimatorSetProxy& e, const World& w) noexcept {
  return transport_method->Transport(*this, e, w);
}

void Particle::Stream(const Real distance) noexcept {
  // update indirect effect due to each perturbation
  for (auto& [perturbation, indirect_effect] : indirect_effects) {
    indirect_effect += perturbation->Stream(*this, distance);
  }
  // update position
  position += direction * distance;
}

void Particle::Scatter(const Real& mu, const Energy& e) noexcept {
  // update indirect effect due to each perturbation
  for (auto& [perturbation, indirect_effect] : indirect_effects) {
    indirect_effect += perturbation->Scatter(*this, mu, e);
  }
  // update direction and energy
  const Real phi = 2 * constants::pi * std::uniform_real_distribution{}(rng);
  direction = Direction{direction, mu, phi};
  energy = e;
}

Real Particle::GetIndirectEffect(
    const Perturbation* perturbation) const noexcept {
  return indirect_effects.at(perturbation);
}

void Particle::SetPerturbations(const PerturbationSet& perturbations) noexcept {
  indirect_effects = perturbations.GetIndirectEffects();
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

void Particle::MoveSecondariesTo(Bank& bank) noexcept {
  bank += std::move(secondaries);
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
