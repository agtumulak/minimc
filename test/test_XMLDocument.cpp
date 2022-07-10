#include "XMLDocument.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_string.hpp"

TEST_CASE("valid XML file passes") {
  REQUIRE_NOTHROW(XMLDocument{"test_XMLDocument_valid.xml"});
}

TEST_CASE("invalid XML file fails") {
  REQUIRE_THROWS_WITH(
      XMLDocument{"test_XMLDocument_invalid.xml"},
      Catch::Matchers::ContainsSubstring(
          "line 4: column 21\n"
          "no declaration found for element 'unexpected_node'"));
}
