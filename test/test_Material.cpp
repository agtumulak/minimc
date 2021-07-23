#include "Material.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

TEST_CASE("nonexistent material name throws exception") {
  XMLDocument doc{"multigroup.xml"};
  REQUIRE_THROWS_WITH(
      Material::FindNode(doc.root, "nonexistent"),
      Catch::Matchers::Equals(
          "Material node \"nonexistent\" not found. "
          "Must be one of: [\"water\", \"hydrogen\", \"oxygen\", ]"));
}
