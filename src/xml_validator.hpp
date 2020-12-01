#pragma once

#include "xercesc/parsers/XercesDOMParser.hpp"
#include "xercesc/sax/ErrorHandler.hpp"
#include "xercesc/sax/SAXParseException.hpp"

#include <filesystem>
#include <memory>
#include <string>

namespace util {
/// @brief Manages XML parsing, validation, and data extraction
/// @details Most methods are adapted from "C++ Cookbook" by Turkanis,
///          Cogswell, Diggins, and Stephens.
class XMLValidator {
public:
  /// @brief Constructs an XMLValidator from a given file path XMLNode
  XMLValidator(const std::filesystem::path& xml_filepath) noexcept;

private:
  // Xerces' wide-character type
  using XMLChString = std::basic_string<XMLCh>;
  // Helper function for constructing and configuring a XercesDOMParser
  static std::unique_ptr<const xercesc::XercesDOMParser> CreateXercesDOMParser(
      const std::filesystem::path& xml_filepath,
      xercesc::ErrorHandler* const error_handler) noexcept;
  // Convert from XMLCh array to std::string
  static std::string toCharString(const XMLCh* const str);
  // RAII wrapper for initializing and terminating the Xerces system. Thus it
  // is deliberately declared before XMLValidator::parser.
  const class XercesInitializer {
  public:
    // Default constructor. Initializes Xerces.
    XercesInitializer() noexcept;
    // Default destructor. Terminates Xerces.
    ~XercesInitializer() noexcept;
    // Disallow copy/move constructor:
    // https://www.stroustrup.com/C++11FAQ.html#default2
    XercesInitializer(const XercesInitializer&) = delete;
    // Disallow copy/move assignment operator:
    // https://www.stroustrup.com/C++11FAQ.html#default2
    XercesInitializer& operator=(const XercesInitializer&) = delete;
  } init;

  // Manages the lifetime of the error handler passed to XMLValidator::parser.
  // Thus it is deliberately declared before XMLValidator::parser.
  class XercesErrorHandler : public xercesc::ErrorHandler {
  public:
    // Handles warnings encountered during XML parsing
    void warning(const xercesc::SAXParseException& exc) noexcept final;
    // Handles recoverable errors encountered during XML parsing
    void error(const xercesc::SAXParseException& exc) noexcept final;
    // Handles non-recovereable errors encountered during XML parsing
    void fatalError(const xercesc::SAXParseException& exc) noexcept final;
    // Reset this handler object after use
    void resetErrors() noexcept final;
    // Write a warning or error message to stderr. Also prints file location
    // information if available.
    void HandleException(const xercesc::SAXParseException& exc) const noexcept;
  } error_handler;

  // Object that parsed the XML document
  std::unique_ptr<const xercesc::XercesDOMParser> parser;
};
} // namespace util
