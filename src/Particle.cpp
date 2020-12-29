#include "Particle.hpp"

// Particle

//// public

Particle::Particle() noexcept
    : position{0., 0., 0.}, direction{1., 0., 0.}, energy{Group{1}} {}

Particle::operator bool() const noexcept { return alive; }

const Point& Particle::GetPosition() const noexcept { return position; };

const Point& Particle::GetDirection() const noexcept { return direction; };

const Energy& Particle::GetEnergy() const noexcept { return energy; };

const Cell& Particle::GetCell() const {
  if (cell) {
    return *cell;
  }
  throw std::runtime_error("Particle does not belong to a Cell");
}

void Particle::SetCell(const Cell& c) noexcept { cell = &c; };
