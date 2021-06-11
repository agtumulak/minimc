#pragma once

#include <filesystem>
#include <vector>

/// @brief Loads an HDF5 file created by pandas
/// @details The pandas DataFrame must be formatted as follows:
///          - There is one and only one column.
///          - The index is a MultiIndex with @f$ N @f$ levels.
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

private:
  // Helper function to read HDF5 file from pandas and return flattened array
  static std::vector<std::vector<double>>
  ReadPandasAxis(const std::filesystem::path& hdf5_filepath);
  // Helper function to read flattened array from HDF5 file
  std::vector<double>
  ReadPandasValues(const std::filesystem::path& hdf5_filepath);
  // Outer vector contains each axis. Inner vector contains values of an axis.
  std::vector<std::vector<double>> axes;
  // The values of the rectangular grid are stored in a flattened array.
  std::vector<double> values;
};
