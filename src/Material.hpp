#pragma once

#include "BasicTypes.hpp"
#include "Nuclide.hpp"
#include "pugixml.hpp"

#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

/// @brief Aggregates nuclear and physical properties
/// @details Nuclear properties are in CrossSection objects
class Material {
public:
  /// @brief Cross sections are Real numbers
  using CrossSection = Real;
  /// @brief Find `material` node with given name in XML document
  /// @param root Root node of existing XML document
  /// @param material_name Value of `name` attribute of `material` node
  /// @exception std::runtime_error `material` node with matching `name`
  ///            attribute not found
  static pugi::xml_node
  FindNode(const pugi::xml_node& root, const std::string& material_name);
  /// @brief Constructs a Material from a material node of an XML document
  /// @param root Root node of existing XML document
  /// @param name Value of `name` attribute of `material` node
  /// @param all_nuclides Nuclide objects which may appear in this Material
  Material(
      const pugi::xml_node& root, const std::string& name,
      const std::vector<std::shared_ptr<const Nuclide>>& all_nuclides);
  /// @brief Sample the distance a Particle will travel before colliding
  Real SampleCollisionDistance(
      std::minstd_rand& rng, const Particle& p) const noexcept;
  /// @brief Sample a Nuclide given that a Particle has collided inside this
  ///        Material
  const Nuclide&
  SampleNuclide(std::minstd_rand& rng, const Particle& p) const noexcept;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;
  /// @brief Atomic number density
  const Real number_density;

private:
  using AtomFractionMap = std::map<std::shared_ptr<const Nuclide>, Real>;
  // Helper function to assign all Nuclide objects which make up Material with
  // the given name. Nuclide objects from all_nuclides are shared between all
  // Material objects.
  static AtomFractionMap AssignNuclides(
      const pugi::xml_node& root, const std::string& name,
      const std::vector<std::shared_ptr<const Nuclide>>& all_nuclides);
  // Return the total <em>microscopic</em> cross section
  CrossSection GetMicroscopicTotal(const Particle& p) const noexcept;
  // All Nuclide objects which make up this Material and their associated atom
  // fractions
  const AtomFractionMap afracs;
};
