#include "World.hpp"

#include <algorithm>

// World

//// public

World::World(const pugi::xml_node& root)
    : surfaces{CreateCSGSurfaces(root)}, nuclides{CreateNuclides(root)},
      materials{CreateMaterials(root, nuclides)}, cells{CreateCells(
                                                      root, surfaces,
                                                      materials)} {}

const Cell& World::FindCellContaining(const Point& p) const {
  const auto cell_it =
      std::find_if(cells.begin(), cells.end(), [&p](const Cell& cell) {
        return cell.Contains(p);
      });
  if (cell_it != cells.end()) {
    return *cell_it;
  }
  throw std::runtime_error("Point does not belong to any Cell");
}

//// private

std::vector<std::shared_ptr<const CSGSurface>>
World::CreateCSGSurfaces(const pugi::xml_node& root) {
  std::vector<std::shared_ptr<const CSGSurface>> all_surfaces;
  for (const auto& cell_node : root.child("cells")) {
    for (const auto& surface_node : cell_node) {
      auto surface_name = surface_node.attribute("name").as_string();
      auto surface_it = std::find_if(
          all_surfaces.cbegin(), all_surfaces.cend(),
          [&surface_name](const std::shared_ptr<const CSGSurface>& surface) {
            return surface->name == surface_name;
          });
      if (surface_it == all_surfaces.cend()) {
        all_surfaces.push_back(CSGSurface::Create(root, surface_name));
      }
    }
  }
  return all_surfaces;
}

std::vector<std::shared_ptr<const Nuclide>>
World::CreateNuclides(const pugi::xml_node& root) {
  std::vector<std::shared_ptr<const Nuclide>> all_nuclides;
  for (const auto& cell_node : root.child("cells")) {
    auto material_name = cell_node.attribute("material").as_string();
    auto material_node = Material::FindNode(root, material_name);
    for (const auto& material_nuclide_node : material_node) {
      auto nuclide_name = material_nuclide_node.attribute("name").as_string();
      auto nuclide_it = std::find_if(
          all_nuclides.cbegin(), all_nuclides.cend(),
          [&nuclide_name](const std::shared_ptr<const Nuclide>& nuclide) {
            return nuclide->name == nuclide_name;
          });
      if (nuclide_it == all_nuclides.cend()) {
        all_nuclides.push_back(
            std::make_shared<const Nuclide>(root, nuclide_name));
      }
    }
  }
  return all_nuclides;
}

std::vector<std::shared_ptr<const Material>> World::CreateMaterials(
    const pugi::xml_node& root,
    const std::vector<std::shared_ptr<const Nuclide>>& all_nuclides) {
  std::vector<std::shared_ptr<const Material>> all_materials;
  for (const auto& cell_node : root.child("cells")) {
    auto material_name = cell_node.attribute("material").as_string();
    auto material_it = std::find_if(
        all_materials.cbegin(), all_materials.cend(),
        [&material_name](const std::shared_ptr<const Material>& material) {
          return material->name == material_name;
        });
    if (material_it == all_materials.cend()) {
      all_materials.push_back(
          std::make_shared<const Material>(root, material_name, all_nuclides));
    }
  }
  return all_materials;
}

std::vector<Cell> World::CreateCells(
    const pugi::xml_node& root,
    const std::vector<std::shared_ptr<const CSGSurface>>& all_surfaces,
    const std::vector<std::shared_ptr<const Material>>&
        all_materials) noexcept {
  std::vector<Cell> cells;
  for (const auto& cell_node : root.child("cells")) {
    cells.emplace_back(
        root, cell_node.attribute("name").as_string(), all_surfaces,
        all_materials);
  }
  return cells;
}
