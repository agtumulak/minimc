#include "Cell.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

// Cell

//// public

Cell::Cell(
    const pugi::xml_node& root, const std::string& cell_name,
    const CSGSurfaceVector& all_surfaces,
    const MaterialVector& all_materials) noexcept
    : name{cell_name}, surface_senses{AssignSurfaceSenses(
                           root, cell_name, all_surfaces)},
      material{AssignMaterial(root, cell_name, all_materials)} {}

bool Cell::Contains(const Point& p) const noexcept {
  return std::all_of(
      surface_senses.cbegin(), surface_senses.cend(),
      [&p](const auto& surface_ptr_bool_pair) {
        const auto& [surface_ptr, is_surface_senses] = surface_ptr_bool_pair;
        // Are the Cell and the Point on the same side of the current surface?
        return surface_ptr->Contains(p) == is_surface_senses;
      });
}

std::tuple<std::shared_ptr<const CSGSurface>, Real>
Cell::NearestSurface(const Point& p, const Point& d) const {
  Real min_distance = std::numeric_limits<Real>::infinity();
  const auto& nearest_it = std::min_element(
      surface_senses.cbegin(), surface_senses.cend(),
      [&p, &d, &min_distance](const auto& lhs, const auto& rhs) {
        Real lhs_distance = lhs.first->Distance(p, d);
        if (lhs_distance < rhs.first->Distance(p, d)) {
          min_distance = lhs_distance;
          return true;
        }
        return false;
      });
  if (nearest_it == surface_senses.cend()) {
    throw std::runtime_error(
        "Particle in Cell \"" + name + "\" could not find nearest surface");
  }
  return std::make_tuple((*nearest_it).first, min_distance);
}

Real Cell::SampleCollisionDistance(
    std::minstd_rand& rng, const Particle& p) const noexcept {
  return material->SampleCollisionDistance(rng, p);
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

std::shared_ptr<const Material> Cell::AssignMaterial(
    const pugi::xml_node& root, const std::string& cell_name,
    const MaterialVector& all_materials) noexcept {
  // void cells do not have a material
  if (cell_name.empty()) {
    return nullptr;
  }
  auto material_name = root.child("cells")
                           .find_child_by_attribute("name", cell_name.c_str())
                           .attribute("material")
                           .as_string();
  auto material_it = std::find_if(
      all_materials.cbegin(), all_materials.cend(),
      [&material_name](const std::shared_ptr<const Material>& material) {
        return material->name == material_name;
      });
  // World::CreateMaterials should have created all Material objects needed
  assert(material_it != all_materials.cend());
  return *material_it;
};
