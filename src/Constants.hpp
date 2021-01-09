#pragma once

#define _USE_MATH_DEFINES

#include "BasicTypes.hpp"

#include <cmath>

namespace constants {
/// @brief Fundamental mathematical constant
/// @details Deprecate with
///          <a href="https://en.cppreference.com/w/cpp/numeric/constants">
///          C++20 Mathematical Constants</a>
constexpr Real pi = M_PI;
} // namespace constants
