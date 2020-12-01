#pragma once

#include <filesystem>

/// @file

/// @brief Validates an XML file.
/// @exception std::runtime_error Parsing failed. Also gives file location and
///            message information if available.
void ValidateXML(const std::filesystem::path& xml_filepath);
