#pragma once

#include "BasicTypes.hpp"
#include "ContinuousMap.hpp"

namespace pugi {
class xml_node;
}

/// @brief Encapsulates all cross section data as well as associated metadata
struct ContinuousEvaluation {
  /// @brief Constructs ContinuousEvaluation from an evaluation node of an XML
  ///        document
  ContinuousEvaluation(const pugi::xml_node& evaluation_node);
  /// @brief Returns true if the relative difference between this evaluation's
  ///        temperature and some requested temperature is sufficiently small
  ///        enough
  bool IsValid(Temperature requested) const noexcept;
  /// @brief Pointwise values of the evaluation
  const ContinuousMap<ContinuousEnergy, MicroscopicCrossSection> xs;
  /// @brief The temperature at which the data was evaluated (C++ Core
  ///        Guidelines C.127)
  const Temperature temperature;
};
