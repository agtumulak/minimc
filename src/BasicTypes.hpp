#pragma once

#include <cstdint>
#include <random>
#include <variant>

/// @file

/// @brief Real number @f$ \mathbb{R} @f$ (C++ Core Guidelines P.1)
using Real = double;
/// @brief Energy in MeV
using ContinuousEnergy = Real;
/// @brief Group number for multigroup calculations. Group 1 corresponds to the
///        highest energy. Larger groups correspond to slower energies.
using Group = std::uint64_t;
/// @brief Since the energy type cannot be known until runtime, this type is
///        used by all functions which deal with energy
using Energy = std::variant<ContinuousEnergy, Group>;
/// @brief Typically used to seed a std::linear_congruential_engine
using History = std::minstd_rand::result_type;
