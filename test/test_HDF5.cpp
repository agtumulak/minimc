#include "HDF5DataSet.hpp"
#include "catch2/catch.hpp"

TEST_CASE("valid HDF5 file passes") {
  REQUIRE_NOTHROW(HDF5DataSet{"test_HDF5_valid.hdf5"});
}

TEST_CASE("invalid HDF5 file fails") {
  REQUIRE_THROWS(HDF5DataSet{"test_HDF5_invalid.hdf5"});
}
