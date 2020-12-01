#include "catch2/catch.hpp"
#include "xml_validator.hpp"

TEST_CASE("parse simple XML file"){
  util::XMLValidator xml_file("test/test_xerces_helpers.xml");
  REQUIRE(true);
}
