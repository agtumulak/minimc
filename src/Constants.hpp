#pragma once

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
/// @brief Used to determine whether a Direction is too close to an axis. A
///        value of 0 means all directions are too close to the axis. A value
///        of 1 means any direction is far away enough from the axis. Must be
///        in (0, 1).
constexpr Real on_axis_tolerance = 0.9;
/// @brief The temperature I personally find most comfortable
constexpr Temperature room_temperature = 293.6;
/// @brief Boltzmann constant in MeV per kelvin
constexpr Real boltzmann = 8.617333262145e-11;
/// @brief Neutron mass in convenient units: MeV * (cm / s)^-2
constexpr Real neutron_mass = 1.045354912280858e-18;
/// @brief Relative temperature difference below which temperatures are to be
///        considered equal for cross section evaluation purposes
constexpr Real relative_temperature_difference_tolerance = 0.01;
/// @brief Number of times to resample thermal scattering alpha when an
///        unphysical value is encountered
constexpr size_t alpha_resample_limit = 10;
} // namespace constants
