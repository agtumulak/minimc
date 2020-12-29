#include "World.hpp"

#include <algorithm>

// World

//// public

World::World(const pugi::xml_node& root)
    : surfaces{CreateCSGSurfaces(root)}, cells{CreateCells(root, surfaces)} {}

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

std::vector<Cell> World::CreateCells(
    const pugi::xml_node& root,
    const std::vector<std::shared_ptr<const CSGSurface>>&
        all_surfaces) noexcept {
  std::vector<Cell> cells;
  for (const auto& cell_node : root.child("cells")) {
    cells.emplace_back(
        root, cell_node.attribute("name").as_string(), all_surfaces);
  }
  return cells;
}
