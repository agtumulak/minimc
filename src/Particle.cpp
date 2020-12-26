#include "Particle.hpp"

// Particle

//// public

Particle::Particle() noexcept
    : position{0., 0., 0.}, direction{1., 0., 0.}, energy{Group{1}} {}

Particle::operator bool() const noexcept { return alive; }
