#pragma once

#include <iosfwd>

/// @brief Supported particle types
enum class Particle {
  neutron,
  photon,
};

/// @brief Helper function to convert from std::string to Particle
Particle ToParticle(const std::string& name);
