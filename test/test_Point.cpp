#include "Point.hpp"
#include "catch2/catch.hpp"

#include <cmath>

TEST_CASE("overloaded Point binary operators") {
  Point p{1, 1, 0};
  Point q{1, 0, 0};
  // add Point objects
  REQUIRE(p + q == Point{2, 1, 0});
  // subtract Point objects
  REQUIRE(p - q == Point{0, 1, 0});
  // compute scalar product of Point and Real
  REQUIRE(p * 2 == Point{2, 2, 0});
  // compute scalar product of Real and Point
  REQUIRE(2 * p == Point{2, 2, 0});
  // compare inequal Point objects
  REQUIRE_FALSE(p == q);
}

TEST_CASE("Point default constructor") {
  Point p{};
  REQUIRE(p == Point{0, 0, 0});
}

TEST_CASE("member functions work"){
  Point u{1, 1, 0};
  u.Normalize();
  REQUIRE(u == Point{1 / std::sqrt(2), 1 / std::sqrt(2), 0});
  // compute vector functions
  Point p{1, 1, 0};
  Point q{1, 0, 0};
  REQUIRE(p.Dot(q) == 1);
  REQUIRE(p.Cross(q) == Point{0, 0, -1});
}

TEST_CASE("overloaded member operators") {
  Point p{1, 2, 3};
  SECTION("assignment by addition") {
    p += Point{4, 5, 6};
    REQUIRE(p == Point{5, 7, 9});
  }
  SECTION("assignment by division") {
    p /= 2;
    REQUIRE(p == Point{0.5, 1, 1.5});
  }
}

TEST_CASE("Direction componentwise constuctor") {
  Direction d{2, 3, 6}; // a Pythagorean quadruple with diagonal length 7
  REQUIRE(d == Point{2./7., 3./7., 6./7.});
}
