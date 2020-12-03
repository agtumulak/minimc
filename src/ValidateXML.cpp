#include "ValidateXML.hpp"

#include "xercesc/parsers/XercesDOMParser.hpp"
#include "xercesc/sax/ErrorHandler.hpp"
#include "xercesc/sax/SAXParseException.hpp"
#include "xercesc/util/PlatformUtils.hpp"
#include "xercesc/util/XMLString.hpp"

#include <stdexcept>
#include <string>

// Convert from XMLCh array to std::string
std::string toCharString(const XMLCh* const str) {
  // unique_ptr is here because transcode gives us ownership of a char array
  std::unique_ptr<const char[]> p{xercesc::XMLString::transcode(str)};
  return std::string{p.get()};
}

// RAII wrapper for initializing and terminating the Xerces system.
class XercesInitializer {
public:
  // Default constructor. Initializes Xerces.
  XercesInitializer() noexcept { xercesc::XMLPlatformUtils::Initialize(); }
  // Default destructor. Terminates Xerces.
  ~XercesInitializer() noexcept { xercesc::XMLPlatformUtils::Terminate(); }
  // Disallow copy/move constructor:
  // https://www.stroustrup.com/C++11FAQ.html#default2
  XercesInitializer(const XercesInitializer&) = delete;
  // Disallow copy/move assignment operator:
  // https://www.stroustrup.com/C++11FAQ.html#default2
  XercesInitializer& operator=(const XercesInitializer&) = delete;
};

/// Error handler registered to xercesc::XercesDOMParser
class XercesErrorHandler : public xercesc::ErrorHandler {
public:
  // Handles warnings encountered during XML parsing
  void warning(const xercesc::SAXParseException& exc) final {
    HandleException(exc);
  }
  // Handles recoverable errors encountered during XML parsing
  void error(const xercesc::SAXParseException& exc) final {
    HandleException(exc);
  }
  // Handles non-recovereable errors encountered during XML parsing
  void fatalError(const xercesc::SAXParseException& exc) final {
    HandleException(exc);
  }
  // Reset this handler object after use
  void resetErrors() noexcept final {}
  // Write a warning or error message to stderr.
  void HandleException(const xercesc::SAXParseException& exc) const {
    std::string error_message;
    auto filename = toCharString(exc.getSystemId());
    auto line_number = exc.getLineNumber();
    auto column_number = exc.getColumnNumber();
    auto have_file_location =
        !filename.empty() && line_number != 0 && column_number != 0;
    if (have_file_location) {
      error_message += filename + ": line " + std::to_string(line_number) +
                       ": column " + std::to_string(column_number) + "\n";
    }
    error_message += toCharString(exc.getMessage()) + "\n";
    throw std::runtime_error(error_message);
  }
};

void ValidateXML(const std::filesystem::path& xml_filepath) {
  XercesInitializer init;
  XercesErrorHandler error_handler;
  xercesc::XercesDOMParser parser;
  parser.setValidationScheme(xercesc::XercesDOMParser::Val_Always);
  parser.setDoSchema(true);
  parser.setDoNamespaces(true);
  parser.setErrorHandler(&error_handler);
  parser.parse(xml_filepath.string().c_str());
}
