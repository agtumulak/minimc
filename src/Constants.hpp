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
} // namespace constants
