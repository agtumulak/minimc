#pragma once

#include <iosfwd>

/// @brief All possible (mutually exclusive, so no total) states that a Particle
///        may undergo
/// @details Assumption of mutual exclusitivity is used when sampling a
///          ContinuousReaction from Continuous::reactions
enum class Reaction {
  birth,
  capture,
  scatter,
  fission,
  leak,
};

/// @brief Helper function to convert from std::string to Reaction
Reaction ToReaction(const std::string& name);
