#include "Point.hpp"
#include "catch2/catch.hpp"

TEST_CASE("overloaded Point binary operators") {
  Point p{1, 1, 0};
  Point q{1, 0, 0};
  // add Point objects
  REQUIRE(p + q == Point{2, 1, 0});
  // subtract Point objects
  REQUIRE(p - q == Point{0, 1, 0});
  // compute inner product of Point objects
  REQUIRE(p * q == 1);
  // compare inequal Point objects
  REQUIRE_FALSE(p == q);
}
