#pragma once

#include "BasicTypes.hpp"

#include <cmath>
#include <cstddef>

/// @brief Returns a value of epsilon for Catch2 Approx assertion macro
/// @details The value of epsilon is such that the mean of i.i.d samples
///          falling within 3 standard deviations (99.7%) of the mean value is
///          accepted.
namespace epsilon {
/// @brief Returns epsilon for Bernoulli distribution
/// @param p Probability of success
/// @param samples Sample size of experiment
inline Real Bernoulli(Real p, size_t samples) {
  if (p == 0){
    return 0;
  }
  return 3 * std::sqrt((1 - p) / (p * samples));
}

/// @brief Returns epsilon for exponential distribution
/// @param samples Epsilon only depends on the number of samples
inline Real Exponential(size_t samples){
  return 3 * std::sqrt(1. / samples);
}

} // namespace epsilon
