#pragma once

#define _USE_MATH_DEFINES

#include "BasicTypes.hpp"

#include <cmath>
#include <limits>

namespace constants {
/// @brief Fundamental mathematical constant
/// @details Deprecate with
///          <a href="https://en.cppreference.com/w/cpp/numeric/constants">
///          C++20 Mathematical Constants</a>
constexpr Real pi = M_PI;
/// @brief Additional distance to stream Particle to cross CSGSurface properly
constexpr Real nudge = 10 * std::numeric_limits<Real>::epsilon();
/// @brief The temperature I personally find most comfortable
constexpr Temperature room_temperature = 293.6;
/// @brief Boltzmann constant in MeV per kelvin
constexpr Real boltzmann = 8.617333262145e-11;
} // namespace constants
