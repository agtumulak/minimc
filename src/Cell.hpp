#pragma once

#include "CSGSurface.hpp"
#include "pugixml.hpp"

#include <map>
#include <memory>
#include <string>

class Particle;

/// @brief A subset of @f$ \mathbb{R}^{3} @f$ defined by constructive solid
///        geometry (CSG) surfaces
class Cell {
private:
  // A polymorphic vector of CSGSurface objects. Multiple Cell objects can
  // share the same CSGSurface object.
  using CSGSurfaceVector = std::vector<std::shared_ptr<const CSGSurface>>;
  // Map from (polymorphic) CSGSurface object to a bool corresponding to the
  // "sense" of the CSGSurface. Multiple Cell objects can share the same
  // CSGSurface object.
  using SurfaceSenses = std::map<std::shared_ptr<const CSGSurface>, bool>;

public:
  /// @brief Constructs a Cell from an XML document
  /// @param root Root node of existing XML document
  /// @param cell_name Value of `name` attribute of `cell` node in XML document
  /// @param all_surfaces CSGSurface objects which may appear in this Cell
  Cell(
      const pugi::xml_node& root, const std::string& cell_name,
      const CSGSurfaceVector& all_surfaces) noexcept;
  /// @brief Returns true if the point lies inside this Cell
  bool Contains(const Point& p) const noexcept;
  /// @brief Returns true if both Cell objects are the same object
  bool operator==(const Cell& rhs) const noexcept;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;
  /// @brief All CSGSurface objects which make up this Cell.
  /// @details All CSGSurfaces have a CSGSurface::Contains method which returns
  ///          `true` if a Point lies on the "negative" side of the CSGSurface.
  ///          This map returns `true` if the Cell lies in the "negative" side
  ///          of the CSGSurface and `false` otherwise.
  const SurfaceSenses surface_senses;

private:
  // Helper function to assign all CSGSurface objects which make up Cell with
  // the given name. CSGSurface objects from all_surfaces are shared between
  // all Cell objects.
  static SurfaceSenses AssignSurfaceSenses(
      const pugi::xml_node& root, const std::string& cell_name,
      const CSGSurfaceVector& all_surfaces) noexcept;
};
