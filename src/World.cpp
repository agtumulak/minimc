#include "World.hpp"

#include "CSGSurface.hpp"
#include "Constants.hpp"
#include "Material.hpp"
#include "Nuclide.hpp"
#include "ScalarField.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>

// World

//// public

World::World(const pugi::xml_node& root)
    : surfaces{CreateCSGSurfaces(root)}, nuclides{CreateNuclides(root)},
      materials{CreateMaterials(root, nuclides)},
      temperature{CreateTemperature(root.child("temperature"))},
      cells{CreateCells(root, surfaces, materials, temperature)} {}

const Cell& World::FindCellContaining(const Point& p) const {
  const auto cell_it =
      std::find_if(cells.begin(), cells.end(), [&p](const Cell& cell) {
        return cell.Contains(p);
      });
  if (cell_it != cells.end()) {
    return *cell_it;
  }
  throw std::runtime_error(
      "Point does not belong to any Cell. Please check"
      "that all space is either assigned a material or void.");
}

bool World::HasConstantTemperature() const noexcept {
  return temperature->IsConstant();
}

std::shared_ptr<const CSGSurface>
World::FindSurfaceByName(const std::string& name) const {
  const auto surface_it = std::find_if(
      surfaces.cbegin(), surfaces.cend(),
      [&name](const auto surface_ptr) { return surface_ptr->name == name; });
  if (surface_it == surfaces.cend()) {
    throw std::runtime_error(
        "Surface \"" + name + "\" not found. Must be one of: [" +
        std::accumulate(
            surfaces.cbegin(), surfaces.cend(), std::string{},
            [](const auto& accumulated, const auto surface_ptr) noexcept {
              return accumulated + "\"" + surface_ptr->name + "\", ";
            }) +
        "]");
  }
  else {
    return *surface_it;
  }
}

std::shared_ptr<const Nuclide>
World::FindNuclideByName(const std::string& name) const {
  const auto nuclide_it = std::find_if(
      nuclides.cbegin(), nuclides.cend(),
      [&name](const auto nuclide_ptr) { return nuclide_ptr->name == name; });
  if (nuclide_it == nuclides.cend()) {
    throw std::runtime_error(
        "Nuclide \"" + name + "\" not found. Must be one of: [" +
        std::accumulate(
            nuclides.cbegin(), nuclides.cend(), std::string{},
            [](const auto& accumulated, const auto nuclide_ptr) noexcept {
              return accumulated + "\"" + nuclide_ptr->name + "\", ";
            }) +
        "]");
  }
  else {
    return *nuclide_it;
  }
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
    if (std::string{cell_node.name()} == "void") {
      // Skip void cells
      continue;
    }
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
        // Find nuclide node and construct Nuclide
        const auto& energy_type_node{root.child("nuclides").first_child()};
        const auto& nuclide_node{
            energy_type_node.find_child_by_attribute("name", nuclide_name)};
        if (!nuclide_node) {
          throw std::runtime_error(
              std::string{"Nuclide node \""} + nuclide_name +
              "\" not found. Must be one of: [" +
              std::accumulate(
                  energy_type_node.begin(), energy_type_node.end(),
                  std::string{},
                  [](const auto& accumulated,
                     const auto& nuclide_node) noexcept {
                    return accumulated + "\"" +
                           nuclide_node.attribute("name").as_string() + "\", ";
                  }) +
              "]");
        }
        all_nuclides.push_back(std::make_shared<const Nuclide>(nuclide_node));
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
    if (std::string{cell_node.name()} == "void") {
      // Skip void cells
      continue;
    }
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

std::unique_ptr<const ScalarField>
World::CreateTemperature(const pugi::xml_node& temperature_node) {
  if (!temperature_node) {
    // default temperature is room temperature
    return std::make_unique<const ConstantField>(constants::room_temperature);
  }
  return ScalarField::Create(temperature_node.first_child());
}

std::vector<Cell> World::CreateCells(
    const pugi::xml_node& root,
    const std::vector<std::shared_ptr<const CSGSurface>>& all_surfaces,
    const std::vector<std::shared_ptr<const Material>>& all_materials,
    const std::shared_ptr<const ScalarField> temperature) noexcept {
  std::vector<Cell> cells;
  for (const auto& cell_node : root.child("cells")) {
    cells.emplace_back(cell_node, all_surfaces, all_materials, temperature);
  }
  return cells;
}
