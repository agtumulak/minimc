#include "Cell.hpp"

#include <algorithm>

// Cell

//// public

Cell::Cell(
    const pugi::xml_node& root, const std::string& cell_name,
    const CSGSurfaceVector& all_surfaces) noexcept
    : name{cell_name}, surface_senses{
                           AssignSurfaceSenses(root, name, all_surfaces)} {}

bool Cell::Contains(const Point& p) const noexcept {
  return std::all_of(
      surface_senses.cbegin(), surface_senses.cend(),
      [&p](const auto& surface_ptr_bool_pair) {
        const auto& [surface_ptr, is_surface_senses] = surface_ptr_bool_pair;
        // Are the Cell and the Point on the same side of the current surface?
        return surface_ptr->Contains(p) == is_surface_senses;
      });
}

bool Cell::operator==(const Cell& rhs) const noexcept { return this == &rhs; }

//// private

Cell::SurfaceSenses Cell::AssignSurfaceSenses(
    const pugi::xml_node& root, const std::string& cell_name,
    const CSGSurfaceVector& all_surfaces) noexcept {
  SurfaceSenses cell_surfaces;
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

    const std::string sense{surface_node.attribute("sense").as_string()};
    bool is_within;
    if (sense == "-1") {
      is_within = true;
    }
    else if (sense == "+1") {
      is_within = false;
    }
    else {
      assert(false); // this should have been caught by the validator
    }
    cell_surfaces.emplace(*surface_it, is_within);
  }
  return cell_surfaces;
}
