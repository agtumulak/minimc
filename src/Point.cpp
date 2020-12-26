#include "Point.hpp"

Point operator+(const Point& lhs, const Point& rhs) noexcept {
  return Point(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
}

Point operator-(const Point& lhs, const Point& rhs) noexcept {
  return Point(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
}

Real operator*(const Point& lhs, const Point& rhs) noexcept {
  return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

bool operator==(const Point& lhs, const Point& rhs) noexcept {
  return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

Point::Point(const pugi::xml_node& pointtype_node) noexcept
    : x{pointtype_node.attribute("x").as_double()},
      y{pointtype_node.attribute("y").as_double()},
      z{pointtype_node.attribute("z").as_double()} {}

Point::Point(const Real& x, const Real& y, const Real& z) noexcept
    : x{x}, y{y}, z{z} {}
