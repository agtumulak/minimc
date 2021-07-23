#include "HDF5DataSet.hpp"

#include "H5Cpp.h"

#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <string>

// HDF5DataSet

//// public

// https://support.hdfgroup.org/HDF5/doc/cpplus_RM/readdata_8cpp-example.html
HDF5DataSet::HDF5DataSet(const std::filesystem::path& hdf5_filepath)
    : axes{ReadPandasAxis(hdf5_filepath)},
      values{ReadPandasValues(hdf5_filepath)}, strides{
                                                   ComputeAxisStrides(axes)} {}

const std::vector<double>& HDF5DataSet::GetAxis(size_t level) const noexcept {
  return axes.at(level);
}

//// private

std::vector<std::vector<double>>
HDF5DataSet::ReadPandasAxis(const std::filesystem::path& hdf5_filepath) {
  std::vector<std::vector<double>> result;
  const H5::H5File file{hdf5_filepath, H5F_ACC_RDONLY};
  // Read number of levels in pandas MultiIndex
  const auto pandas_group = file.openGroup("/pandas");
  size_t levels; pandas_group.openAttribute("axis1_nlevels")
      .read(H5::PredType::NATIVE_ULONG, &levels);
  // Read each axis into result
  for (size_t i = 0; i < levels; i++) {
    const auto level_dataset = pandas_group.openDataSet(
        std::string{"axis1_level"} + std::to_string(i));
    result.emplace_back(
        static_cast<size_t>(level_dataset.getSpace().getSimpleExtentNpoints()),
        0);
    level_dataset.read(result.back().data(), H5::PredType::NATIVE_DOUBLE);
  }
  return result;
}

std::vector<double>
HDF5DataSet::ReadPandasValues(const std::filesystem::path& hdf5_filepath) {
  std::vector<double> result;
  const H5::H5File file{hdf5_filepath, H5F_ACC_RDONLY};
  // Read elements into result
  const auto block_dataset =
      file.openGroup("/pandas").openDataSet("block0_values");
  result = std::vector<double>(
      static_cast<size_t>(block_dataset.getSpace().getSimpleExtentNpoints()),
      0);
  block_dataset.read(result.data(), H5::PredType::NATIVE_DOUBLE);
  return result;
}

std::vector<size_t> HDF5DataSet::ComputeAxisStrides(
    const std::vector<std::vector<double>>& axes) noexcept {
  // Initialize result. `axes` is assumed to have size one or greater.
  std::vector<size_t> result(axes.size());
  // Element `i` in `result` contains the size of axis `i+1`
  std::transform(
      std::next(axes.cbegin(), 1), axes.cend(), result.begin(),
      [](const auto& axis) { return axis.size(); });
  // Last element is 1 since the stride at the innermost axis (level `N`) will
  // be one.
  result.back() = 1;
  // Element `i` in `result` contains the product of the sizes of axes `i+1`,
  // `i+1`, ..., `N`.
  std::partial_sum(
      result.crbegin(), result.crend(), result.rbegin(), std::multiplies());
  return result;
}
