#include "HDF5DataSet.hpp"
#include "catch2/catch.hpp"

TEST_CASE("valid HDF5 file passes") {
  REQUIRE_NOTHROW(HDF5DataSet{"test_HDF5_valid.hdf5"});
}

TEST_CASE("invalid HDF5 file fails") {
  REQUIRE_THROWS(HDF5DataSet{"test_HDF5_invalid.hdf5"});
}

TEST_CASE("element access using variadic member function works"){
  HDF5DataSet hdf5_dataset{"test_HDF5_valid.hdf5"};
  // std::upper_bound is used so these will return the value at the next-highest
  // data point, (3.81, 18.92, 573.6)
  REQUIRE(hdf5_dataset.at(0.00, 0.0005, 293.6) == Approx(1.142070e-2));
  REQUIRE(hdf5_dataset.at(3.8099, 18.9199, 573.599) == Approx(1.142070e-2));
  // std::upper_bound is used so this will return the highest data point
  REQUIRE(hdf5_dataset.at(35.30, 18.92, 293.6) == Approx(7.921300e-17));
  REQUIRE(hdf5_dataset.at(93.99, 310.99, 573.599) == Approx(7.921300e-17));
  // maximum value of zeroth axis is 94.0
  REQUIRE_THROWS_WITH(
      hdf5_dataset.at(94.0, 0.0005, 293.6), "Upper bound not found");
  // maximum value of last axis is 573.6
  REQUIRE_THROWS_WITH(
      hdf5_dataset.at(0.00, 0.0005, 573.6), "Upper bound not found");
}
