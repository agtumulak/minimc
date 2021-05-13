#pragma once

#include <string>

/// @brief All possible (mutually-exclusive, so no total) reactions
///        regardless of incident particle type
enum class Reaction {
  capture,
  scatter,
  fission,
};

/// @brief Helper function to convert from std::string to Reaction
Reaction ToReaction(const std::string& name) noexcept;
