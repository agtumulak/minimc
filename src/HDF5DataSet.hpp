#pragma once

#include "BasicTypes.hpp"

#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

/// @brief Loads an HDF5 file created by pandas
/// @details The pandas DataFrame must be formatted as follows:
///          - There is one and only one column.
///          - The index is a MultiIndex with @f$ N @f$ levels. It must be a
///            MultiIndex even if has only one level.
///          - Each level of the MultiIndex corresponds to one axis of a
///            multidimensional array.
///          - Each level of the MultiIndex is sorted.
///          - The first level of the MultiIndex (level @f$ 0 @f$) changes the slowest.
///          - The last level of the MultiIndex (level @f$ N - 1 @f$) changes
///            the fastest.
///          - All data types are of `double` type.
///          The desired ordering of index levels can be achieved with
///          <a href="https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.sort_index.html">
///          `pandas.DataFrame.sort_index(level=[`axis0`, ..., `axisN`])`</a>.
///          - <tt>key='pandas'</tt> when
///            <a href="https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.to_hdf.html">
///            `pandas.DataFrame.to_hdf`</a> is called.
class HDF5DataSet {
public:
  /// @brief Constructs an HDF5DataSet from an hdf5 file
  HDF5DataSet(const std::filesystem::path& hdf5_filepath);
  /// @brief Returns the value at the given axis indices
  /// @details The user is responsible for knowing what each level of the HDF5
  ///          DataFrame corresponds to. For example, consider a DataFrame of
  ///          space-dependent temperatures with level 0 corresponding to
  ///          x-position, level 1 corresponding to y-position, and level 2
  ///          corresponding to z-position. One can retrieve the temperature at
  ///          index 4 of level 0, index 2 of level 1, and index 0 of level 2
  ///          with `at(4, 2, 0)`.
  template <typename... Args> Real at(size_t index, Args... inner_indices) {
    return at_from_base(0, 0, index, inner_indices...);
  }

private:
  // Helper function to read HDF5 file from pandas and return flattened array
  static std::vector<std::vector<double>>
  ReadPandasAxis(const std::filesystem::path& hdf5_filepath);
  // Helper function to read flattened array from HDF5 file
  static std::vector<double>
  ReadPandasValues(const std::filesystem::path& hdf5_filepath);
  // Helper function to precompute strides
  static std::vector<size_t>
  ComputeAxisStrides(const std::vector<std::vector<double>>& axes) noexcept;
  // Starting from a given base offset, using innermost levels starting at
  // given level, return value from full flattened array
  template <typename... Args>
  Real
  at_from_base(size_t base, size_t level, size_t index, Args... inner_indices) {
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
  std::vector<std::vector<double>> axes;
  // The values of the rectangular grid are stored in a flattened array.
  std::vector<double> values;
  // Distance of a single step in each level
  std::vector<size_t> strides;
};
