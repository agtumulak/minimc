#include "xml_validator.hpp"

#include "xercesc/util/PlatformUtils.hpp"
#include "xercesc/util/XMLString.hpp"

#include <iostream>

using util::XMLValidator;

// XercesErrorHandler

void XMLValidator::XercesErrorHandler::warning(
    const xercesc::SAXParseException& exc) noexcept {
  std::cerr << "Warning: ";
  HandleException(exc);
};

void XMLValidator::XercesErrorHandler::error(
    const xercesc::SAXParseException& exc) noexcept {
  std::cerr << "Error: ";
  HandleException(exc);
};

void XMLValidator::XercesErrorHandler::fatalError(
    const xercesc::SAXParseException& exc) noexcept {
  std::cerr << "Fatal Error: ";
  HandleException(exc);
};

void XMLValidator::XercesErrorHandler::resetErrors() noexcept {};

void XMLValidator::XercesErrorHandler::HandleException(
    const xercesc::SAXParseException& exc) const noexcept {
  auto filename = toCharString(exc.getSystemId());
  auto line_number = exc.getLineNumber();
  auto column_number = exc.getColumnNumber();
  auto have_file_location =
      !filename.empty() && line_number != 0 && column_number != 0;
  if (have_file_location) {
    std::cerr << filename << ": line " << line_number << ": column "
              << column_number << std::endl;
  }
  std::cerr << toCharString(exc.getMessage()) << std::endl;
}

// XMLValidator

XMLValidator::XMLValidator(const std::filesystem::path& xml_filepath) noexcept
    : parser(CreateXercesDOMParser(xml_filepath, &error_handler)) {}

std::unique_ptr<const xercesc::XercesDOMParser>
XMLValidator::CreateXercesDOMParser(
    const std::filesystem::path& xml_filepath,
    xercesc::ErrorHandler* const error_handler) noexcept {
  auto parser = std::make_unique<xercesc::XercesDOMParser>();
  parser->setValidationScheme(xercesc::XercesDOMParser::Val_Always);
  parser->setDoSchema(true);
  parser->setDoNamespaces(true);
  parser->setErrorHandler(error_handler);
  parser->parse(xml_filepath.c_str());
  return parser;
};

std::string XMLValidator::toCharString(const XMLCh* const str) {
  // unique_ptr is here because transcode gives us ownership of a char array
  std::unique_ptr<const char[]> p{xercesc::XMLString::transcode(str)};
  return std::string{p.get()};
}

// XercesInitializer

XMLValidator::XercesInitializer::XercesInitializer() noexcept {
  xercesc::XMLPlatformUtils::Initialize();
}

XMLValidator::XercesInitializer::~XercesInitializer() noexcept {
  xercesc::XMLPlatformUtils::Terminate();
}
