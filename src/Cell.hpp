#pragma once

#include "CSGSurface.hpp"
#include "Material.hpp"
#include "pugixml.hpp"

#include <map>
#include <memory>
#include <random>
#include <string>
#include <tuple>

class Particle;

/// @brief A subset of @f$ \mathbb{R}^{3} @f$ defined by constructive solid
///        geometry (CSG) surfaces
class Cell {
private:
  // A polymorphic vector of CSGSurface objects. Multiple Cell objects can
  // share the same CSGSurface object.
  using CSGSurfaceVector = std::vector<std::shared_ptr<const CSGSurface>>;
  // Map from (polymorphic) CSGSurface object to a bool corresponding to the
  // "sense" of the CSGSurface. Multiple Cell objects can share the same
  // CSGSurface object.
  using SurfaceSenses = std::map<std::shared_ptr<const CSGSurface>, bool>;
  // A vector of Material objects. Multiple Cell objects can share the same
  // Material object.
  using MaterialVector = std::vector<std::shared_ptr<const Material>>;

public:
  /// @brief Constructs a Cell from an XML document
  /// @param root Root node of existing XML document
  /// @param cell_name Value of `name` attribute of `cell` node in XML document
  /// @param all_surfaces CSGSurface objects which may appear in this Cell
  /// @param all_materials Material objects which may appear in this Cell
  Cell(
      const pugi::xml_node& root, const std::string& cell_name,
      const CSGSurfaceVector& all_surfaces,
      const MaterialVector& all_materials) noexcept;
  /// @brief Returns true if the point lies inside this Cell
  bool Contains(const Point& p) const noexcept;
  /// @brief Finds the nearest CSGSurface from a point along a given direction
  /// @param p position in the Cell
  /// @param d direction to search along
  /// @exception std::runtime_error No CSGSurface found
  /// @return Returns a tuple to nearest CSGSurface and distance to it
  std::tuple<std::shared_ptr<const CSGSurface>, Real>
  NearestSurface(const Point& p, const Direction& d) const;
  /// @brief Sample the distance a particle will travel before colliding
  Real SampleCollisionDistance(RNG& rng, const Particle& p) const noexcept;
  /// @brief Returns true if both Cell objects are the same object
  bool operator==(const Cell& rhs) const noexcept;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;
  /// @brief All CSGSurface objects which make up this Cell.
  /// @details All CSGSurfaces have a CSGSurface::Contains method which returns
  ///          `true` if a Point lies on the "negative" side of the CSGSurface.
  ///          This map returns `true` if the Cell lies in the "negative" side
  ///          of the CSGSurface and `false` otherwise.
  const SurfaceSenses surface_senses;
  /// @brief Material this Cell is made of. A nullptr corresponds to a void.
  const std::shared_ptr<const Material> material;

private:
  // Helper function to assign all CSGSurface objects which make up Cell with
  // the given name. CSGSurface objects from all_surfaces are shared between
  // all Cell objects.
  static SurfaceSenses AssignSurfaceSenses(
      const pugi::xml_node& root, const std::string& cell_name,
      const CSGSurfaceVector& all_surfaces) noexcept;
  // Helper function to assign Material object used by the Cell with the given
  // name. Material objects from all_materials may be shared between Cell
  // objects.
  static std::shared_ptr<const Material> AssignMaterial(
      const pugi::xml_node& root, const std::string& cell_name,
      const MaterialVector& all_materials) noexcept;
};
