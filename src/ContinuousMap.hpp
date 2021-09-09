#pragma once

#include "BasicTypes.hpp"

#include <map>
#include <utility>

/// @brief Continuously maps elements from a domain to a range, provided a
///        limited set of points, by interpolating
/// @tparam Key Type of the domain
/// @tparam T Type of the range
template <typename Key, typename T> class ContinuousMap {
public:
  /// @brief Type used to store elements internally
  using elements_type = std::map<Key, T>;

  /// @brief Constructs ContinuousMap by assigning elements directly
  ContinuousMap(elements_type&& other) : elements{std::move(other)} {}

  /// @brief Returns a const reference to the value at a given key
  const T& at(const Key k) const noexcept {
    // TODO: Interpolate and handle edge cases
    return elements.upper_bound(k)->second;
  }

protected:
  /// @brief This class essentially wraps an STL container
  const elements_type elements;
};

/// @brief Like Map, but stores elements as the CDF of some random variable.
///        Stores CDF values as keys so that std::map::upper_bound() can be
///        used.
/// @tparam T The type of the random variable
template <typename T> class CDF : public ContinuousMap<Real, T> {
public:
  /// @brief Constructs a CDF from a std::map
  CDF(typename ContinuousMap<Real, T>::elements_type&& other)
      : ContinuousMap<Real, T>{std::move(other)} {};

  /// @brief Samples a value from the CDF and returns the sampled key
  const T& Sample(RNG& rng) const noexcept {
    return this->elements.upper_bound(std::uniform_real_distribution{}(rng))
        ->second;
  }
};
