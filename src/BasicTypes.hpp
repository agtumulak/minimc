#pragma once

#include <cstdint>
#include <random>
#include <variant>

/// @file

/// @brief Index for Bins classes
using BinIndex = size_t;
/// @brief Real number @f$ \mathbb{R} @f$ (C++ Core Guidelines P.1)
using Real = double;
/// @brief See @ref estimators_scoring_functions
using Score = Real;
/// @brief Temperature in kelvins
using Temperature = Real;
/// @brief Microscopic cross section in barns
using MicroscopicCrossSection = Real;
/// @brief Macroscopic cross section in inverse cm
using MacroscopicCrossSection = Real;
/// @brief Energy in MeV
using ContinuousEnergy = Real;
/// @brief Group number for multigroup calculations. Group 1 corresponds to the
///        highest energy. Larger groups correspond to slower energies.
using Group = std::uint64_t;
/// @brief Since the energy type cannot be known until runtime, this type is
///        used by all functions which deal with energy
using Energy = std::variant<ContinuousEnergy, Group>;
/// @brief Users may want to use a different random number generator in the
///        future
using RNG = std::minstd_rand;
