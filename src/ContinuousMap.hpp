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
  /// @brief Default constructs an empty ContinuousMap
  ContinuousMap() {}
  /// @brief Constructs ContinuousMap by assigning elements directly
  ContinuousMap(elements_type&& other) : elements{std::move(other)} {}
  /// @brief Returns a const reference to the value at a given key
  /// @todo Interpolate and handle edge cases
  const T& at(const Key k) const noexcept {
    return elements.upper_bound(k)->second;
  }
  /// @brief Returns a reference to the value at a given key
  T& operator[](const Key& k) { return elements[k]; }

protected:
  /// @brief This class essentially wraps an STL container
  elements_type elements;
};

/// @brief Contains a nested ContinuousMap type
/// @details https://stackoverflow.com/a/1500289/5101335
/// @tparam D Nest depth
template <size_t D> struct NestedContinuousMap {
  /// @brief Type contained by a ContinuousMap when there are D levels of
  ///        nesting
  using MappedType =
      ContinuousMap<Real, typename NestedContinuousMap<D - 1>::MappedType>;
};

/// @brief Contains the base case for NestedContinuousMap
template <> struct NestedContinuousMap<0> {
  /// @brief Type contained by a ContinuousMap when there is no nesting
  using MappedType = Real;
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
