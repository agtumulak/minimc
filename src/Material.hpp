#pragma once

#include "BasicTypes.hpp"
#include "pugixml.hpp"

#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <vector>

class Nuclide;
class Particle;

/// @brief Aggregates nuclear and physical properties
/// @details Nuclear properties are in NuclearData objects
class Material {
public:
  /// @brief Nuclide objects and their associated atom fractions
  using AtomFractions = std::map<std::shared_ptr<const Nuclide>, Real>;
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
  /// @brief Return the majorant <em>microscopic</em> cross section
  /// @details This is not to be confused with the <em>global</em> majorant
  ///          across all Materials in the problem
  MicroscopicCrossSection
  GetMicroscopicMajorant(const Particle& p) const noexcept;
  /// @brief Return the total <em>microscopic</em> cross section
  MicroscopicCrossSection GetMicroscopicTotal(const Particle& p) const noexcept;
  /// @brief Unique, user-defined identifier (C++ Core Guidelines C.131)
  const std::string name;
  /// @brief All Nuclide objects which make up this Material and their
  ///        associated atom fractions (C++ Core Guidelines C.131)
  const AtomFractions afracs;
  /// @brief Atomic number density (C++ Core Guidelines C.131)
  const Real number_density;

private:
  // Helper function to assign all Nuclide objects which make up Material with
  // the given name. Nuclide objects from all_nuclides are shared between all
  // Material objects.
  static AtomFractions AssignNuclides(
      const pugi::xml_node& root, const std::string& name,
      const std::vector<std::shared_ptr<const Nuclide>>& all_nuclides);

};
