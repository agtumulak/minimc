#pragma once

#include "CSGSurface.hpp"
#include "pugixml.hpp"

#include <memory>
#include <string>
#include <vector>

/// @brief A subset of @f$ \mathbb{R}^{3} @f$ defined by constructive solid
///        geometry (CSG) surfaces
class Cell {
private:
  // A polymorphic vector of CSGSurface objects. Multiple Cell objects can
  // share the same CSGSurface object.
  using CSGSurfaceVector = std::vector<std::shared_ptr<const CSGSurface>>;

public:
  /// @brief Constructs a Cell from an XML document
  /// @param root Root node of existing XML document
  /// @param cell_name Value of `name` attribute of `cell` node in XML document
  /// @param all_surfaces CSGSurface objects which may appear in this Cell
  Cell(
      const pugi::xml_node& root, const std::string& cell_name,
      const CSGSurfaceVector& all_surfaces) noexcept;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;

private:
  // Helper function to assign all CSGSurface objects which make up Cell with
  // the given name. CSGSurface objects from all_surfaces are shared between
  // all Cell objects.
  static CSGSurfaceVector AssignCSGSurfaces(
      const pugi::xml_node& root, const std::string& cell_name,
      const CSGSurfaceVector& all_surfaces) noexcept;
  // All CSGSurface objects which make up this Cell
  const CSGSurfaceVector surfaces;
};
