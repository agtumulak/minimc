#pragma once

#include "BasicTypes.hpp"
#include "Point.hpp"
#include "pugixml.hpp"

#include <memory>
#include <string>

/// @brief A surface to be used in constructive solid geometry
class CSGSurface {
public:
  /// @brief Factory method to create new CSGSurface from an XML document
  /// @param root Root node of existing XML document
  /// @param name Name of the CSGSurface in the XML document
  /// @returns A `std::unique_ptr` to the constructed CSGSurface (C++ Core
  ///          Guidelines R.30)
  /// @exception std::runtime_error `surface` node with matching `name`
  ///            attribute not found
  /// @todo Throw exception when there exist more than one `surface` node with
  ///       the same `name` attribute
  static std::unique_ptr<const CSGSurface>
  Create(const pugi::xml_node& root, const std::string& name);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~CSGSurface() noexcept;
  /// @brief Returns true if Point is in the "negative" side of the (open)
  ///        surface. A point <em>on</em> the surface or outside is considered
  ///        to be on the "positive" side.
  virtual bool Contains(const Point& p) const noexcept = 0;
  /// @brief Return the distance from a given origin Point to the CSGSurface
  ///        along a given Direction
  /// @returns A positive finite distance if Direction is pointed towards the
  ///          CSGSurface, otherwise positive infinity
  virtual Real
  Distance(const Point& origin, const Direction& direction) const noexcept = 0;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;

protected:
  /// @brief Specialized form of quadratic equation solver
  /// @details Used for CSGSurface intersection distance calculations. Returns
  ///          least positive real solution if it exists; infinity otherwise.
  static Real SolveQuadratic(Real a, Real b, Real c) noexcept;
  /// @brief Constructs a CSGSurface from a surface node of an XML document
  /// @param surface_node The requested `surface` node in the XML document
  CSGSurface(const pugi::xml_node& surface_node) noexcept;

private:
  // Possible derived CSGSurface types
  enum class SurfaceType {
    Sphere,
    PlaneX,
    CylinderX,
  };
  // Converts std::string to CSGSurface::SurfaceType
  static SurfaceType ToSurfaceType(const std::string& surface_name) noexcept;
};

/// @brief A sphere. What else has to be said?
class Sphere : public CSGSurface {
public:
  /// @brief Constructs a Sphere from a `sphere` node of an XML document
  /// @param sphere_node The requested `sphere` node in the XML document
  Sphere(const pugi::xml_node& sphere_node) noexcept;
  /// @brief Returns the distance from a given origin Point to the Sphere along
  ///        a given Direction using the
  ///        <a href="https://en.wikipedia.org/wiki/Line-sphere_intersection">
  ///        line-sphere intersection algorithm</a>
  /// @param origin Starting point from where distance will be calculated.
  /// @param direction Unit vector. Must be normalized to unity.
  Real Distance(
      const Point& origin, const Direction& direction) const noexcept override;
  /// @brief Implements CSGSurface method
  bool Contains(const Point& p) const noexcept override;

private:
  const Point center;
  const Real radius;
};

/// @brief Plane perpendicular to the x-axis
class PlaneX : public CSGSurface {
public:
  /// @brief Constructs a PlaneX from a `planex` node of an XML document
  /// @param planex_node The requested `planex` node in the XML document
  PlaneX(const pugi::xml_node& planex_node) noexcept;

  /// Returns the distance from a given origin Point to the PlaneX along a
  /// given Direction
  Real Distance(
      const Point& origin, const Direction& direction) const noexcept override;
  /// @brief Implements CSGSurface method
  bool Contains(const Point& p) const noexcept override;

private:
  // completely defines the plane x - c = 0
  const Real c;
};

/// @brief Infinitely long cylinder coaxial with the x-axis
class CylinderX : public CSGSurface {
public:
  /// @brief Constructs a CylinderX from a `cylinderx` node of an XML document
  /// @param cylinderx_node The requested `cylinderx` node in the XML document
  CylinderX(const pugi::xml_node& cylinderx_node) noexcept;
  /// @brief Returns the distance from a given origin Point to the CylinderX
  ///        along a given Direction
  Real Distance(
      const Point& origin, const Direction& direction) const noexcept override;
  /// @brief IMplements CSGSurface method
  bool Contains(const Point& p) const noexcept override;

private:
  // completely defines the cylinder r^2 = y^2 + z^2
  const Real radius;
};
