#include "Particle.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "Material.hpp"
#include "Nuclide.hpp"
#include "Reaction.hpp"
#include "World.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <random>
#include <stdexcept>

// Particle::TransportOutcome

//// public

Particle::TransportOutcome&
Particle::TransportOutcome::operator+=(const TransportOutcome& rhs) noexcept {
  estimator += rhs.estimator;
  std::move(
      rhs.banked.cbegin(), rhs.banked.cend(),
      std::back_insert_iterator(banked));
  return *this;
}

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

Particle::TransportOutcome Particle::Transport(const World& w) noexcept {
  TransportOutcome result;
  SetCell(w.FindCellContaining(position));
  while (alive) {
    const auto collision = SampleCollisionDistance();
    const auto [nearest_surface, surface_crossing] =
        cell->NearestSurface(position, direction);
    if (collision < surface_crossing) {
      result.estimator.at(Estimator::Event::collision) += 1;
      Stream(collision);
      const auto& nuclide = SampleNuclide();
      result.estimator.at(Estimator::Event::implicit_fission) +=
          nuclide.GetNuBar(*this) *
          nuclide.GetReaction(*this, Reaction::fission) /
          nuclide.GetTotal(*this);
      nuclide.Interact(*this);
    }
    else {
      result.estimator.at(Estimator::Event::surface_crossing) += 1;
      Stream(surface_crossing + constants::nudge);
      SetCell(w.FindCellContaining(position));
      if (!cell->material) {
        Kill();
      }
    }
  }
  return result;
}

void Particle::Kill() noexcept { alive = false; }

void Particle::Stream(const Real distance) noexcept {
  position += direction * distance;
}

const Point& Particle::GetPosition() const noexcept { return position; };

void Particle::SetDirectionIsotropic() noexcept {
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

Real Particle::SampleCollisionDistance() noexcept {
  return std::exponential_distribution{
      cell->material->number_density *
      cell->material->GetMicroscopicTotal(*this)}(rng);
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
