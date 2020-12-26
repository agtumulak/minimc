#include "Cell.hpp"

#include <algorithm>

// Cell

//// public

Cell::Cell(
    const pugi::xml_node& root, const std::string& name,
    const CSGSurfaceVector& all_surfaces) noexcept
    : surfaces{AssignCSGSurfaces(root, name, all_surfaces)}, name{name} {}

//// private

Cell::CSGSurfaceVector Cell::AssignCSGSurfaces(
    const pugi::xml_node& root, const std::string& cell_name,
    const CSGSurfaceVector& all_surfaces) noexcept {
  CSGSurfaceVector cell_surfaces;
  for (const auto& surface_node :
       root.child("cells").find_child_by_attribute("name", cell_name.c_str())) {
    const auto surface_name = surface_node.attribute("name").as_string();
    const auto surface_it = std::find_if(
        all_surfaces.cbegin(), all_surfaces.cend(),
        [&surface_name](const std::shared_ptr<const CSGSurface>& surface) {
          return surface->name == surface_name;
        });
    // World::CreateCSGSurfaces should have created all CSGSurface objects
    // needed
    assert(surface_it != all_surfaces.cend());
    cell_surfaces.push_back(*surface_it);
  }
  return cell_surfaces;
}
