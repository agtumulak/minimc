#pragma once

#include "CSGSurface.hpp"
#include "Cell.hpp"
#include "pugixml.hpp"

#include <memory>
#include <vector>

/// @brief Represents the state of a nuclear system
/// @details Encapsualtes geometry, materials, and sources
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
  // Helper function to create vector of cells
  static std::vector<Cell> CreateCells(
      const pugi::xml_node& root,
      const std::vector<std::shared_ptr<const CSGSurface>>&
          all_surfaces) noexcept;

  // All unique (polymorphic) CSGSurface objects shared among Cell objects
  const std::vector<std::shared_ptr<const CSGSurface>> surfaces;

public:
  /// @brief Cells that appear in the World (C++ Core Guidelines C.131)
  const std::vector<Cell> cells;
};
