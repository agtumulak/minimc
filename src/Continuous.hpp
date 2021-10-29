#pragma once

#include "BasicTypes.hpp"
#include "ContinuousEvaluation.hpp"
#include "ContinuousReaction.hpp"
#include "Interaction.hpp"
#include "pugixml.hpp"

#include <memory>
#include <vector>

class Particle;

/// @brief Contains cross sections which are indexed by continuous energy values
class Continuous : public Interaction {
public:
  /// @brief Constructs continuous energy nuclear data from a particle node of
  ///        an XML document
  Continuous(const pugi::xml_node& particle_node);
  /// @brief Returns the largest cross section that may be found within the
  ///        current Cell
  /// @details Currently this is the cross section at the majorant temperature
  ///          in the Cell
  MicroscopicCrossSection
  GetMajorant(const Particle& p) const noexcept override;
  /// @brief Returns the total cross section for a given Particle
  /// @details This is not guaranteed to be consistent with the sum of all
  ///          mutually exclusive reactions. The total cross section is meant
  ///          to be a user-provided quantity to speed up calculations.
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept override;
  /// @brief Interact with a Particle, updating its state
  void Interact(Particle& p) const noexcept override;

private:
  // Helper function for reaction cross section construction
  static std::vector<std::unique_ptr<const ContinuousReaction>>
  CreateReactions(const pugi::xml_node& particle_node);
  // Returns true if any reaction modifies the total cross section even if the
  // total cross section was evaluated at the target temperature
  bool ReactionsModifyTotal(const Particle& p) const noexcept;
  // Cross section data for mutually exclusive reactions
  const std::vector<std::unique_ptr<const ContinuousReaction>> reactions;
  // Total cross section provided in nuclear data files
  const ContinuousEvaluation total;
};
