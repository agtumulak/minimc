#include "BasicTypes.hpp"
#include "CSGSurface.hpp"
#include "Point.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

#include <limits>
#include <memory>

TEST_CASE("nonexistent surface name throws exception") {
  XMLDocument doc{"simple_multigroup.xml"};
  REQUIRE_THROWS_WITH(
      CSGSurface::Create(doc.root, "nonexistent"),
      "Surface node \"nonexistent\" not found. Must be one of: "
      "[\"inner sphere\", \"middle sphere\", \"outer sphere\", ]");
}

TEST_CASE("compute distances to Sphere") {
  XMLDocument doc{"simple_multigroup.xml"};
  const std::unique_ptr<const CSGSurface> surface{
      CSGSurface::Create(doc.root, "inner sphere")};

  // discriminant less than zero; no intersection
  REQUIRE(
      surface->Distance(Point{2, 0, 0}, Direction{0, 0, 1}) ==
      std::numeric_limits<Real>::infinity());
  // discriminant equal zero; no intersection; grazing
  REQUIRE(
      surface->Distance(Point{1, 0, 0}, Direction{0, 0, 1}) ==
      std::numeric_limits<Real>::infinity());
  // outside sphere; headed towards sphere
  REQUIRE(surface->Distance(Point{0, 0, -2}, Direction{0, 0, 1}) == Approx(1));
  // on sphere; leaving sphere
  REQUIRE(surface->Distance(Point{0, 0, -1}, Direction{0, 0, 1}) == Approx(2));
  // in sphere; leaving sphere
  REQUIRE(surface->Distance(Point{0, 0, 0}, Direction{0, 0, 1}) == Approx(1));
  // on sphere; headed away from sphere
  REQUIRE(
      surface->Distance(Point{0, 0, 1}, Direction{0, 0, 1}) ==
      std::numeric_limits<Real>::infinity());
  // outside sphere; headed away from sphere
  REQUIRE(
      surface->Distance(Point{0, 0, 2}, Direction{0, 0, 1}) ==
      std::numeric_limits<Real>::infinity());
}
