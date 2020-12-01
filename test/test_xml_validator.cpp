#include "catch2/catch.hpp"
#include "xml_validator.hpp"

TEST_CASE("valid XML file passes") {
  REQUIRE_NOTHROW(ValidateXML("test/test_xml_validator_valid.xml"));
}

TEST_CASE("invalid XML file fails") {
  REQUIRE_THROWS_WITH(
      ValidateXML("test/test_xml_validator_invalid.xml"),
      Catch::Matchers::Contains(
          "line 4: column 21\n"
          "no declaration found for element 'unexpected_node'"));
}
