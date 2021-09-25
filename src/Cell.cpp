#include "Cell.hpp"

#include "CSGSurface.hpp"
#include "Material.hpp"
#include "ScalarField.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string>

// Cell

//// public

Cell::Cell(
    const pugi::xml_node& cell_node, const CSGSurfaceVector& all_surfaces,
    const MaterialVector& all_materials,
    const std::shared_ptr<const ScalarField> global_temperature) noexcept
    : name{cell_node.attribute("name").as_string()},
      surface_senses{AssignSurfaceSenses(cell_node, all_surfaces)},
      material{AssignMaterial(cell_node, all_materials)},
      temperature{AssignTemperature(cell_node, global_temperature)} {}

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
Cell::NearestSurface(const Point& p, const Direction& d) const {
  const auto& nearest_it = std::min_element(
      surface_senses.cbegin(), surface_senses.cend(),
      [&p, &d](const auto& lhs, const auto& rhs) {
        return lhs.first->Distance(p, d) < rhs.first->Distance(p, d);
      });
  if (nearest_it == surface_senses.cend()) {
    throw std::runtime_error(
        "Particle in Cell \"" + name + "\" could not find nearest surface");
  }
  const auto& nearest_surface = (*nearest_it).first;
  // TODO: Cache result of nearest surface distance instead of computing again
  return std::make_tuple(nearest_surface, nearest_surface->Distance(p, d));
}

bool Cell::operator==(const Cell& rhs) const noexcept { return this == &rhs; }

//// private

Cell::SurfaceSenses Cell::AssignSurfaceSenses(
    const pugi::xml_node& cell_node,
    const CSGSurfaceVector& all_surfaces) noexcept {
  SurfaceSenses cell_surfaces;
  for (const auto& surface_node : cell_node) {
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
    const pugi::xml_node& cell_node,
    const MaterialVector& all_materials) noexcept {
  // void cells do not have a material
  if (cell_node.attribute("name").empty()) {
    return nullptr;
  }
  auto material_name = cell_node.attribute("material").as_string();
  auto material_it = std::find_if(
      all_materials.cbegin(), all_materials.cend(),
      [&material_name](const std::shared_ptr<const Material>& material) {
        return material->name == material_name;
      });
  // World::CreateMaterials should have created all Material objects needed
  assert(material_it != all_materials.cend());
  return *material_it;
};

const std::shared_ptr<const ScalarField> Cell::AssignTemperature(
    const pugi::xml_node& cell_node,
    const std::shared_ptr<const ScalarField> global_temperature) noexcept {
  return cell_node.attribute("temperature")
             ? std::make_shared<const ConstantField>(
                   cell_node.attribute("temperature").as_double())
             : global_temperature;
}
