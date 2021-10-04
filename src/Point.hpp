#pragma once

#include "BasicTypes.hpp"
#include "pugixml.hpp"

/// @brief Point in @f$ \mathbb{R}^{3} @f$
class Point {
  /// @brief Returns the vector sum of two Point objects
  friend Point operator+(const Point& lhs, const Point& rhs) noexcept;
  /// @brief Returns the vector difference of two Point objects
  friend Point operator-(const Point& lhs, const Point& rhs) noexcept;
  /// @brief Returns Point with each element multiplied by rhs
  friend Point operator*(const Point& lhs, const Real& rhs) noexcept;
  /// @brief Returns Point with each element multiplied by lhs
  friend Point operator*(const Real& lhs, const Point& rhs) noexcept;
  /// @brief Returns Point with each element divided by rhs
  friend Point operator/(const Point& lhs, const Real& rhs) noexcept;

public:
  /// @brief Default constructor. Creates a Point at the origin;
  Point() noexcept;
  /// @brief Constructs a Point from a `PointType` node
  /// @details `PointType` is a `complexType` defined in the minimc XML schema.
  Point(const pugi::xml_node& pointtype_node) noexcept;
  /// @brief Constructs a Point with the given components
  Point(const Real& x, const Real& y, const Real& z) noexcept;
  /// @brief Scales the point to satisfy @f$ \lVert v \rVert = 1 @f$
  const Point& Normalize() noexcept;
  /// @brief Returns the inner product of this and another Point
  Real Dot(const Point& rhs) const noexcept;
  /// @brief Returns the cross product of this and another Point
  Point Cross(const Point& rhs) const noexcept;
  /// @brief Adds the given Point to the current Point
  Point& operator+=(const Point& rhs) noexcept;
  /// @brief Divides each element of this Point with rhs
  Point& operator/=(const Real& rhs) noexcept;
  /// @brief Returns true if each corresponding element is equal
  bool operator==(const Point& rhs) const noexcept;

protected:
  /// @brief x component
  Real x{0};
  /// @brief y component
  Real y{0};
  /// @brief z component
  Real z{0};
};

/// @brief Point in @f$ \mathbb{R}^{3} @f$ subject to
///        @f$ \lVert v \rVert = 1 @f$
/// @details All users of this class can assume that the base Point is
///          normalized.
/// @note    Maintainers of this class must take care to keep the underlying
///          Point normalized by calling Normalize() when appropriate.
class Direction : public Point {
public:
  /// @brief Constructs an isotropically sampled Direction
  static Direction CreateIsotropic(RNG& rng) noexcept;
  /// @brief Constructs a new Direction that has a given cosine @f$ \mu @f$
  ///        with respect to an exsiting Direction.
  /// @details Given a Direction @f$ \boldsymbol{\Omega}^{\prime} @f$, the new
  ///          Direction @f$ \boldsymbol{\Omega} @f$ will satisfy @f$
  ///          \boldsymbol{\Omega}^{\prime} \cdot \boldsymbol{\Omega} = \mu
  ///          @f$.
  /// @param d The reference direction @f$ \boldsymbol{\Omega}^{\prime} @f$
  /// @param mu The cosine @f$ \mu @f$.
  /// @param phi The azimuthal angle about @f$ \boldsymbol{\Omega}^{\prime}
  ///            @f$. This is largely arbitrary since scattering is usually
  ///            azimuthally symmetric.
  static Direction CreateAboutDirection(
      const Direction& d, const Real& mu, const Real& phi) noexcept;
  /// @brief Constructs a Direction from a `PointType` node
  /// @details `PointType` is a `complexType` defined in the minimc XML schema.
  Direction(const pugi::xml_node& pointtype_node) noexcept;
  /// @brief Constructs a Direction with the given components
  Direction(const Real& x, const Real& y, const Real& z) noexcept;
  /// @brief Move constructor. Constructs the Direction with the contents of a
  ///        Point. The given Point is normalized.
  Direction(Point&& other) noexcept;
};
