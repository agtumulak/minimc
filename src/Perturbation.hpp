#pragma once

#include <iosfwd>

/// @brief All possible perturbation types
enum class Perturbation {
  total,
  capture,
  scatter,
};

/// @brief Helper function to convert from std::string to Reaction
Perturbation ToPerturbation(const std::string& name);
