#pragma once

#include "CSGSurface.hpp"
#include "Cell.hpp"
#include "Material.hpp"
#include "pugixml.hpp"

#include <memory>
#include <vector>

/// @brief Represents the state of a nuclear system
/// @details Manages the construction of all Surface, Nuclide, Material, and
///          Cell objects. Access to shared members is given by World when
///          necessary.
class World {
public:
  /// @brief Constructs a World from an XML document
  World(const pugi::xml_node& root);
  /// @brief Returns the current Cell occupied by a given Point
  /// @details TODO: More than one Cell may contain a Particle. This should be
  ///          checked during input parsing if possible.
  const Cell& FindCellContaining(const Point& p) const;

private:
  // Helper function to create all CSGSurface objects that appear in Cell
  // objects
  static std::vector<std::shared_ptr<const CSGSurface>>
  CreateCSGSurfaces(const pugi::xml_node& root);
  // Helper function to create all Nuclide objects that appear in all Material
  // objects that appear in Cell objects
  static std::vector<std::shared_ptr<const Nuclide>>
  CreateNuclides(const pugi::xml_node& root);
  // Helper function to create all Material objects that appear in Cell ojects
  static std::vector<std::shared_ptr<const Material>> CreateMaterials(
      const pugi::xml_node& root,
      const std::vector<std::shared_ptr<const Nuclide>>& all_nuclides);
  // Helper function to create vector of cells
  static std::vector<Cell> CreateCells(
      const pugi::xml_node& root,
      const std::vector<std::shared_ptr<const CSGSurface>>& all_surfaces,
      const std::vector<std::shared_ptr<const Material>>&
          all_materials) noexcept;

  // All unique (polymorphic) CSGSurface objects shared among Cell objects
  const std::vector<std::shared_ptr<const CSGSurface>> surfaces;
  // All unique Nuclide objects shared among Material objects
  const std::vector<std::shared_ptr<const Nuclide>> nuclides;
  // All unique Material objects shared among Cell objects
  const std::vector<std::shared_ptr<const Material>> materials;

public:
  /// @brief Cells that appear in the World (C++ Core Guidelines C.131)
  const std::vector<Cell> cells;
};
