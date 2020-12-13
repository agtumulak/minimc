#pragma once

#include "pugixml.hpp"

#include <filesystem>

/// @brief Validates document with Xerces-C++. Loads document with pugixml.
class XMLDocument {
public:
  /// @brief Loads a valid pugixml XML document.
  /// @exception std::runtime_error Validation or DOM loading failed. Gives
  ///            file location and message information if validation failed.
  XMLDocument(const std::filesystem::path& xml_filepath);
  /// @brief Allows access to the loaded pugi::xml_document
  /// @details pugi::xml_node objects are only pointers to the full document.
  ///          Therefore calling code must not let the owning XMLDocument to
  ///          go out of scope.
  pugi::xml_node root;

private:
  // Validates an XML file and returns a path to the XML file.
  static void ValidateXML(const std::filesystem::path& xml_filepath);
  // Owns the entire XML document structure
  pugi::xml_document doc;
};
