#pragma once

#include <iosfwd>

/// @brief All possible (mutually exclusive, so no total) reactions regardless
///        of incident particle type
/// @details Assumption of mutual exclusitivity is used when sampling a
///          ContinuousReaction from Continuous::reactions
enum class Reaction {
  capture,
  scatter,
  fission,
};

/// @brief Helper function to convert from std::string to Reaction
Reaction ToReaction(const std::string& name);
