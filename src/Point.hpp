#pragma once

#include "BasicTypes.hpp"
#include "pugixml.hpp"

/// @brief Point in @f$ \mathbb{R}^{3} @f$
class Point {
  /// @brief Returns the vector sum of two Point objects
  friend Point operator+(const Point& lhs, const Point& rhs) noexcept;
  /// @brief Returns the vector difference of two Point objects
  friend Point operator-(const Point& lhs, const Point& rhs) noexcept;
  /// @brief Returns the inner product of two Point objects
  friend Real operator*(const Point& lhs, const Point& rhs) noexcept;

public:
  /// @brief Constructs a Point from a `PointType` node
  /// @details `PointType` is a `complexType` defined in the minimc XML schema.
  Point(const pugi::xml_node& pointtype_node) noexcept;
  /// @brief Constructs a Point with the given components
  Point(const Real& x, const Real& y, const Real& z) noexcept;
  /// @brief Returns true if each corresponding element is equal
  bool operator==(const Point& rhs) const noexcept;

private:
  /// @brief x component
  Real x{0};
  /// @brief y component
  Real y{0};
  /// @brief z component
  Real z{0};
};
