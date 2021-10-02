#pragma once

#include "BasicTypes.hpp"

#include <iterator>
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
  /// @todo Deprecate in favor of construction from
  ///       HDF5DataSet::ToContinuousMap
  ContinuousMap(elements_type&& other) : elements{std::move(other)} {}
  /// @brief Returns a linearly interpolated value at the given key
  /// @tparam Args Type of inner keys
  /// @todo Support other interpolation methods
  template <typename... Args>
  decltype(auto) at(const Key k, Args... inner_keys) const noexcept {
    const auto it_hi = elements.upper_bound(k);
    if constexpr (sizeof...(Args) == 0){
      // base case
      if (it_hi == elements.cend()){
        // snap to last entry
        return std::prev(it_hi, 1)->second;
      }
      else if (it_hi == elements.cbegin()) {
        // snap to first entry
        return it_hi->second;
      }
      else {
        // linearly interpolate
        const auto& [k_hi, v_hi] = *it_hi;
        const auto& [k_lo, v_lo] = *std::prev(it_hi, 1);
        return v_lo + (v_hi - v_lo) / (k_hi - k_lo) * (k - k_lo);
      }
    }
    else {
      // nested case
      if (it_hi == elements.cend()){
        return std::prev(it_hi, 1)->second.at(inner_keys...);
      }
      else if (it_hi == elements.cbegin()){
        return it_hi->second.at(inner_keys...);
      }
      else {
        const auto& k_hi = it_hi->first;
        const auto& v_hi = it_hi->second.at(inner_keys...);
        const auto& k_lo = std::prev(it_hi, 1)->first;
        const auto& v_lo = std::prev(it_hi, 1)->second.at(inner_keys...);
        return v_lo + (v_hi - v_lo) / (k_hi - k_lo) * (k - k_lo);
      }
    }
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
