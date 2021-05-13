#include "Material.hpp"
#include "Statistics.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

#include <random>

TEST_CASE("nonexistent material name throws exception") {
  XMLDocument doc{"simple_multigroup.xml"};
  REQUIRE_THROWS_WITH(
      Material::FindNode(doc.root, "nonexistent"),
      Catch::Matchers::Equals(
          "Material node \"nonexistent\" not found. "
          "Must be one of: [\"water\", \"hydrogen\", \"oxygen\", ]"));
}
