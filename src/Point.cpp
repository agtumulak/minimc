#include "Point.hpp"

#include "Constants.hpp"

#include <cmath>
#include <random>

// Point

Point operator+(const Point& lhs, const Point& rhs) noexcept {
  return Point{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Point operator-(const Point& lhs, const Point& rhs) noexcept {
  return Point{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

Real operator*(const Point& lhs, const Point& rhs) noexcept {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Point operator*(const Point& lhs, const Real& rhs) noexcept {
  return Point{lhs.x * rhs, lhs.y * rhs, lhs.z * rhs};
}

Point operator*(const Real& lhs, const Point& rhs) noexcept {
  return rhs * lhs;
}

//// public

Point::Point() noexcept {}

Point::Point(const pugi::xml_node& pointtype_node) noexcept
    : x{pointtype_node.attribute("x").as_double()},
      y{pointtype_node.attribute("y").as_double()},
      z{pointtype_node.attribute("z").as_double()} {}

Point::Point(const Real& x, const Real& y, const Real& z) noexcept
    : x{x}, y{y}, z{z} {}

void Point::Normalize() noexcept {
  *this /= std::sqrt(*this * *this); // sorry about this
}

Point& Point::operator+=(const Point& rhs) noexcept {
  x += rhs.x;
  y += rhs.y;
  z += rhs.z;
  return *this;
}

Point& Point::operator/=(const Real& rhs) noexcept {
  x /= rhs;
  y /= rhs;
  z /= rhs;
  return *this;
}

bool Point::operator==(const Point& rhs) const noexcept {
  return x == rhs.x && y == rhs.y && z == rhs.z;
}

// Direction

//// public

Direction Direction::CreateIsotropic(RNG& rng) noexcept {
  Real x = std::uniform_real_distribution{-1., +1.}(rng);
  const Real sin_theta = std::sqrt(1 - x * x);
  const Real phi = std::uniform_real_distribution{0., 2 * constants::pi}(rng);
  Real y = sin_theta * std::cos(phi);
  Real z = sin_theta * std::sin(phi);
  return Direction{x, y, z};
}

Direction::Direction(const pugi::xml_node& pointtype_node) noexcept
    : Point{pointtype_node} {
  Normalize();
}

Direction::Direction(const Real& x, const Real& y, const Real& z) noexcept
    : Point{x, y, z} {
  Normalize();
}
