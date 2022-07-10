#include "Material.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_exception.hpp"

#include <stdexcept>

TEST_CASE("nonexistent material name throws exception") {
  XMLDocument doc{"multigroup.xml"};
  REQUIRE_THROWS_MATCHES(
      Material::FindNode(doc.root, "nonexistent"), std::runtime_error,
      Catch::Matchers::Message(
          "Material node \"nonexistent\" not found. "
          "Must be one of: [\"water\", \"hydrogen\", \"oxygen\", ]"));
}
