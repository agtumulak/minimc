#include "Point.hpp"
#include "catch2/catch.hpp"

TEST_CASE("Point default constructor") {
  Point p{};
  REQUIRE(p == Point{0, 0, 0});
}

TEST_CASE("overloaded Point binary operators") {
  Point p{1, 1, 0};
  Point q{1, 0, 0};
  // add Point objects
  REQUIRE(p + q == Point{2, 1, 0});
  // subtract Point objects
  REQUIRE(p - q == Point{0, 1, 0});
  // compute inner product of Point objects
  REQUIRE(p * q == 1);
  // compute scalar product of Point and Real
  REQUIRE(p * 2 == Point{2, 2, 0});
  // compute scalar product of Real and Point
  REQUIRE(2 * p == Point{2, 2, 0});
  // compare inequal Point objects
  REQUIRE_FALSE(p == q);
}

TEST_CASE("overloaded member operators") {
  Point p{1, 2, 3};
  p += Point{4, 5, 6};
  REQUIRE(p == Point{5, 7, 9});
}
