#pragma once

#include "pugixml.hpp"

#include <filesystem>

/// @brief Validates document with Xerces-C++. Loads document with pugixml.
/// @todo Replace xs:normalizedString by xs:token or more restrictive type in
///       schema
class XMLDocument {
public:
  /// @brief Loads a valid pugixml XML document.
  /// @exception std::runtime_error Validation or DOM loading failed. Gives
  ///            file location and message information if validation failed.
  XMLDocument(const std::filesystem::path& xml_filepath);

private:
  // Validates an XML file with Xerces-C++ and returns a pugi::xml_document
  static pugi::xml_document
  ValidateXML(const std::filesystem::path& xml_filepath);
  // Owns the entire XML document structure
  const pugi::xml_document doc;

public:
  /// @brief Allows access to the root node of XMLDocument::doc
  /// @details pugi::xml_node objects are only pointers to the full document.
  ///          Therefore calling code must not let the owning XMLDocument to
  ///          go out of scope. This member is declared after XMLDocument::doc
  ///          to guarantee construction order.
  const pugi::xml_node root{doc.child("minimc")};
};
