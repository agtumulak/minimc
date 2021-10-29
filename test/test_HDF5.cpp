#include "HDF5DataSet.hpp"
#include "catch2/catch.hpp"

TEST_CASE("valid HDF5 file passes") {
  REQUIRE_NOTHROW(HDF5DataSet<3>{"test_HDF5_valid.hdf5"});
}

TEST_CASE("invalid HDF5 file fails") {
  REQUIRE_THROWS(HDF5DataSet<1>{"test_HDF5_invalid.hdf5"});
}

TEST_CASE("element access using variadic member function works"){
  HDF5DataSet<3> hdf5_dataset{"test_HDF5_valid.hdf5"};
  REQUIRE(hdf5_dataset.at(0, 0, 0) == Approx(4.628770e-02));
  REQUIRE(hdf5_dataset.at(0, 0, 1) == Approx(1.761820e-01));
  REQUIRE(hdf5_dataset.at(0, 1, 0) == Approx(9.757560e-03));
  REQUIRE(hdf5_dataset.at(0, 1, 1) == Approx(2.031040e-02));
  REQUIRE(hdf5_dataset.at(1, 0, 0) == Approx(2.147320e-06));
  REQUIRE(hdf5_dataset.at(1, 0, 1) == Approx(1.212070e-05));
  REQUIRE_THROWS_WITH(
      hdf5_dataset.at(4, 0, 0), "Out-of-range (level, index): (0, 4)");
  REQUIRE_THROWS_WITH(
      hdf5_dataset.at(0, 3, 0), "Out-of-range (level, index): (1, 3)");
  REQUIRE_THROWS_WITH(
      hdf5_dataset.at(0, 0, 2), "Out-of-range (level, index): (2, 2)");
}
