#pragma once

#include "BasicTypes.hpp"
#include "ContinuousMap.hpp"
#include "H5Cpp.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

/// @brief Loads an HDF5 file created by pandas
/// @details The pandas DataFrame must be formatted as follows:
///          - There is one and only one column.
///          - The index is a MultiIndex with @f$ D @f$ levels. It must be a
///            MultiIndex even if has only one level.
///          - Each level of the MultiIndex corresponds to one axis of a
///            multidimensional array.
///          - Each level of the MultiIndex is sorted.
///          - The first level of the MultiIndex (level @f$ 0 @f$) changes the slowest.
///          - The last level of the MultiIndex (level @f$ D - 1 @f$) changes
///            the fastest.
///          - All data types are of `double` type.
///          The desired ordering of index levels can be achieved with
///          <a href="https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.sort_index.html">
///          `pandas.DataFrame.sort_index(level=[`axis0`, ..., `axisD`])`</a>.
///          - <tt>key='pandas'</tt> when
///            <a href="https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_hdf.html">
///            `pandas.DataFrame.to_hdf`</a> is called.
/// @tparam D number of dimensions (independent axes) in HDF5 file
template <size_t D> class HDF5DataSet {
public:
  /// @brief Constructs an HDF5DataSet from an hdf5 file
  /// @details https://support.hdfgroup.org/HDF5/doc/cpplus_RM/readdata_8cpp-example.html
  HDF5DataSet(const std::filesystem::path& hdf5_filepath)
      : axes{ReadPandasAxis(hdf5_filepath)},
        values{ReadPandasValues(hdf5_filepath)}, strides{ComputeAxisStrides(
                                                     axes)} {}
  /// @brief Returns an instance of ContinuousMap with D levels of nesting
  /// @tparam Args Type of outer indices
  /// @param outer_indices Indices of outer levels; populated during recursion
  /// @todo Create ContinuousMap constructor which wraps this function
  template <typename... Args>
  typename NestedContinuousMap<D - sizeof...(Args)>::MappedType
  ToContinuousMap(Args... outer_indices) {
    typename NestedContinuousMap<D - sizeof...(Args)>::MappedType result;
    if constexpr (D - sizeof...(Args) == 0) {
      result = at(outer_indices...);
    }
    else {
      const auto& keys = GetAxis(sizeof...(outer_indices));
      for (size_t index = 0; index < keys.size(); index++) {
        result[keys[index]] = ToContinuousMap(outer_indices..., index);
      }
    }
    return result;
  };
  /// @brief Returns the value at the given axis indices
  /// @details The user is responsible for knowing what each level of the HDF5
  ///          DataFrame corresponds to. For example, consider a DataFrame of
  ///          space-dependent temperatures with level 0 corresponding to
  ///          x-position, level 1 corresponding to y-position, and level 2
  ///          corresponding to z-position. One can retrieve the temperature at
  ///          index 4 of level 0, index 2 of level 1, and index 0 of level 2
  ///          with `at(4, 2, 0)`.
  template <typename... Args>
  Real at(size_t index, Args... inner_indices) const {
    // compile time error if wrong number of indices
    static_assert(1 + sizeof...(inner_indices) == D);
    return at_from_base(0, 0, index, inner_indices...);
  }
  /// @brief Get values of a given axis
  const std::vector<double>& GetAxis(size_t level) const noexcept {
    return axes.at(level);
  }
  /// @brief Returns total number of elements in dataset
  size_t size() const noexcept { return values.size(); }

private:
  // Helper function to read HDF5 file from pandas and return flattened array.
  // Of all the other helper function this is the one that checks for file
  // existence just because it happens to be called first.
  static std::array<std::vector<double>, D>
  ReadPandasAxis(const std::filesystem::path& hdf5_filepath) {
    if (!std::filesystem::exists(hdf5_filepath)) {
      throw std::runtime_error(
          std::string{"File not found: "} + hdf5_filepath.string());
    }
    std::array<std::vector<double>, D> result;
    const H5::H5File file(hdf5_filepath.string(), H5F_ACC_RDONLY);
    // Read number of levels in pandas MultiIndex
    const auto pandas_group = file.openGroup("/pandas");
    size_t levels;
    pandas_group.openAttribute("axis1_nlevels")
        .read(H5::PredType::NATIVE_ULLONG, &levels);
    if (levels != D) {
      throw std::runtime_error(
          hdf5_filepath.string() + ": Expected " + std::to_string(D) +
          " dimensions, but got " + std::to_string(levels));
    }
    // Read each axis into result
    for (size_t i = 0; i < D; i++) {
      const auto level_dataset = pandas_group.openDataSet(
          std::string{"axis1_level"} + std::to_string(i));
      result[i] = std::vector<double>(
          static_cast<size_t>(
              level_dataset.getSpace().getSimpleExtentNpoints()),
          0);
      level_dataset.read(result[i].data(), H5::PredType::NATIVE_DOUBLE);
    }
    return result;
  }
  // Helper function to read flattened array from HDF5 file
  static std::vector<double>
  ReadPandasValues(const std::filesystem::path& hdf5_filepath) {
    std::vector<double> result;
    const H5::H5File file(hdf5_filepath.string(), H5F_ACC_RDONLY);
    // Read elements into result
    const auto block_dataset =
        file.openGroup("/pandas").openDataSet("block0_values");
    result = std::vector<double>(
        static_cast<size_t>(block_dataset.getSpace().getSimpleExtentNpoints()),
        0);
    block_dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);
    return result;
  }
  // Helper function to precompute strides.
  static std::array<size_t, D>
  ComputeAxisStrides(const std::array<std::vector<double>, D>& axes) noexcept {
    // Initialize result
    std::array<size_t, D> result;
    // Element `i` in `result` contains the size of axis `i+1`
    std::transform(
        std::next(axes.cbegin(), 1), axes.cend(), result.begin(),
        [](const auto& axis) { return axis.size(); });
    // Last element is 1 since the stride at the innermost axis (level `D - 1`)
    // will be one.
    result[D - 1] = 1;
    // Compute strides by partially multiplying each level size beginning at
    // level `D - 1`
    std::partial_sum(
        result.crbegin(), result.crend(), result.rbegin(), std::multiplies());
    return result;
  }
  // Starting from a given base offset, using innermost levels starting at
  // given level, return value from full flattened array
  template <typename... Args>
  Real at_from_base(
      size_t base, size_t level, size_t index, Args... inner_indices) const {
    // runtime error if index is out of range
    if (index >= axes.at(level).size()) {
      throw std::out_of_range(
          "Out-of-range (level, index): (" + std::to_string(level) + ", " +
          std::to_string(index) + ")");
    }
    if constexpr (sizeof...(inner_indices) == 0) {
      return values.at(base + index);
    }
    else {
      return at_from_base(
          base + index * strides.at(level), level + 1, inner_indices...);
    }
  }
  // Outer vector contains each axis. Inner vector contains values of an axis.
  const std::array<std::vector<double>, D> axes;
  // The values of the rectangular grid are stored in a flattened array.
  const std::vector<double> values;
  // Distance of a single step in each level
  const std::array<size_t, D> strides;
};
