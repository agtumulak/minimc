#include "BasicTypes.hpp"
#include "CSGSurface.hpp"
#include "Point.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_exception.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <limits>
#include <memory>
#include <stdexcept>

TEST_CASE("nonexistent surface name throws exception") {
  XMLDocument doc{"multigroup.xml"};
  REQUIRE_THROWS_MATCHES(
      CSGSurface::Create(doc.root, "nonexistent"), std::runtime_error,
      Catch::Matchers::Message(
          "Surface node \"nonexistent\" not found. Must be one of: "
          "[\"inner shell\", \"middle shell\", \"outer shell\", ]"));
}

TEST_CASE("compute distances to Sphere") {
  XMLDocument doc{"multigroup.xml"};
  const std::unique_ptr<const CSGSurface> surface{
      CSGSurface::Create(doc.root, "inner shell")};

  // discriminant less than zero; no intersection
  REQUIRE(
      surface->Distance(Point{2, 0, 0}, Direction{0, 0, 1}) ==
      std::numeric_limits<Real>::infinity());
  // discriminant equal zero; no intersection; grazing
  REQUIRE(
      surface->Distance(Point{1, 0, 0}, Direction{0, 0, 1}) ==
      std::numeric_limits<Real>::infinity());
  // outside sphere; headed towards sphere
  REQUIRE_THAT(
      surface->Distance(Point{0, 0, -2}, Direction{0, 0, 1}),
      Catch::Matchers::WithinRel(1.));
  // on sphere; leaving sphere
  REQUIRE_THAT(
      surface->Distance(Point{0, 0, -1}, Direction{0, 0, 1}),
      Catch::Matchers::WithinRel(2.));
  // in sphere; leaving sphere
  REQUIRE_THAT(
      surface->Distance(Point{0, 0, 0}, Direction{0, 0, 1}),
      Catch::Matchers::WithinRel(1.));
  // on sphere; headed away from sphere
  REQUIRE(
      surface->Distance(Point{0, 0, 1}, Direction{0, 0, 1}) ==
      std::numeric_limits<Real>::infinity());
  // outside sphere; headed away from sphere
  REQUIRE(
      surface->Distance(Point{0, 0, 2}, Direction{0, 0, 1}) ==
      std::numeric_limits<Real>::infinity());
}
