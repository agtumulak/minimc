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

Point operator*(const Point& lhs, const Real& rhs) noexcept {
  return Point{lhs.x * rhs, lhs.y * rhs, lhs.z * rhs};
}

Point operator*(const Real& lhs, const Point& rhs) noexcept {
  return rhs * lhs;
}

Point operator/(const Point& lhs, const Real& rhs) noexcept {
  return Point{lhs.x / rhs, lhs.y / rhs, lhs.z / rhs};
}

//// public

Point::Point() noexcept {}

Point::Point(const pugi::xml_node& pointtype_node) noexcept
    : x{pointtype_node.attribute("x").as_double()},
      y{pointtype_node.attribute("y").as_double()},
      z{pointtype_node.attribute("z").as_double()} {}

Point::Point(const Real& x, const Real& y, const Real& z) noexcept
    : x{x}, y{y}, z{z} {}

const Point& Point::Normalize() noexcept {
  *this /= std::sqrt(Dot(*this));
  return *this;
}

Real Point::Dot(const Point& rhs) const noexcept {
  return x * rhs.x + y * rhs.y + z * rhs.z;
}

Point Point::Cross(const Point& rhs) const noexcept {
  return Point{
      y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x};
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

Direction::Direction(
    const Direction& d, const Real& mu, const Real& phi) noexcept {
  // Determine if d is outside an "hourglass" of directions which we consider
  // "too close" to the x axis. Note edge case where on_axis_tolerance = 1 will
  // still reject d.x = - 1 but on_axis_tolerance should never be set to 1
  // anyways.
  const bool is_off_xaxis = d.x <= constants::on_axis_tolerance &&
                           d.x > -constants::on_axis_tolerance;
  // Construct another Point `x` which is not parallel to `d`. This will be a
  // unit vector along x (or y).
  const auto xaxis = is_off_xaxis ? Point{1, 0, 0} : Point{0, 1, 0};
  // Take cross product, the returned unit vector is in the y-z (or x-z)
  // plane and orthogonal to `d`
  const Direction u{d.Cross(xaxis)};
  // Take cross product, the returned vector is also orthogonal to `u` and `d`
  const Direction v{d.Cross(u)};
  // `u` and `v` span a plane of points which are orthogonal to `d`. Use the
  // three vectors to construct the new direction
  const auto u_comp = std::sqrt(1 - mu * mu) * std::cos(phi) * u;
  const auto v_comp = std::sqrt(1 - mu * mu) * std::sin(phi) * v;
  const auto d_comp = mu * d;
  const Direction omega{u_comp + v_comp + d_comp};
  x = omega.x;
  y = omega.y;
  z = omega.z;
  // No need to Normalize()
}

Direction::Direction(Point&& other) noexcept : Point{std::move(other)} {
  Normalize();
}
