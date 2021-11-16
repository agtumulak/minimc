#pragma once

#include "BasicTypes.hpp"
#include "Interaction.hpp"
#include "Reaction.hpp"

#include <map>
#include <optional>
#include <vector>

namespace pugi {
class xml_node;
}
struct State;

/// @brief Contains cross sections which are indexed by discrete energy groups
/// @details Groups are integers in `[1,G]`. Group `1` corresponds to the
///          highest energy. Group `G` corresponds to the lowest energy.
/// @todo Refactor Multigroup::reactions into vector instead of map as is now
///       done with Continuous::reactions
class Multigroup : public Interaction {
public:
  /// @brief Constructs multigroup nuclear data from a particle node of an XML
  ///        document
  /// @details Creates cross sections for each Reaction. If no corresponding
  ///          reaction node is found, it is not added to the map.
  /// @param particle_node The requested particle node in the XML document
  /// @exception std::runtime_error Number of entries is not consistent with
  ///            number of groups
  Multigroup(const pugi::xml_node& particle_node);
  /// @brief Returns the largest cross section that may be found within the
  ///        current Cell for a given particle State
  /// @details Currently this is the same as GetTotal().
  MicroscopicCrossSection GetMajorant(const State& s) const noexcept override;
  /// @brief Returns the total cross section for a given particle State
  MicroscopicCrossSection GetTotal(const State& s) const noexcept override;
  /// @brief Interact with the particle, updating its State
  void Interact(State& s) const noexcept override;

private:
  // Groupwise cross sections indexed by one Group
  class OneDimensional {
  private:
    // The STL container underlying a OneDimensional object
    using elements_type = std::vector<Real>;

  public:
    // Constructs Multigroup::OneDimensional from vector of Real
    OneDimensional(const elements_type& elements);
    // Constructs Multigroup::OneDimensional from a `GroupXS` type node. Refer
    // to XML schema for structure of `GroupXS` type node.
    OneDimensional(const pugi::xml_node& groupxs_node);
    // Returns a const reference to the `g`-th group
    const Real& at(const Group g) const;
    // Returns an iterator corresponding to most energetic Group
    elements_type::iterator begin() noexcept;
    // Returns a const_iterator corresponding to most energetic Group
    elements_type::const_iterator begin() const noexcept;
    // Returns an iterator corresponding to the Group following the least
    // energetic Group
    elements_type::iterator end() noexcept;
    // Returns a const_iterator corresponding to the Group following the least
    // energetic Group
    elements_type::const_iterator end() const noexcept;

  private:
    elements_type elements;
  };
  // Groupwise cross sections indexed by two Groups
  class TwoDimensional {
  protected:
    // The STL container underlying a TwoDimensional object
    using elements_type = std::vector<OneDimensional>;

  public:
    // Constructs Multigroup::TwoDimensional from a `GroupXS` type node. Refer
    // to XML schema for structure of `GroupXS` type node.
    TwoDimensional(const pugi::xml_node& groupxs_node);
    // Returns a const reference to the `g`-th group
    const OneDimensional& at(const Group g) const;
    // Returns an iterator corresponding to most energetic Group
    elements_type::const_iterator begin() const noexcept;
    // Returns an iterator corresponding to the Group following the least
    // energetic Group
    elements_type::const_iterator end() const noexcept;

  protected:
    elements_type elements;
  };
  // Like TwoDimensional, but each OneDimensional element is normalized
  class NormalizedTwoDimensional : public TwoDimensional {
  public:
    // Constructs Multigroup::TwoDimensional from a `GroupXS` type node. Refer
    // to XML schema for structure of `GroupXS` type node.
    NormalizedTwoDimensional(const pugi::xml_node& groupxs_node);
  };
  using ReactionsMap = std::map<Reaction, OneDimensional>;
  // Returns the user-specified number of groups under the `multigroup` node
  static Group GroupStructureSize(const pugi::xml_node& root) noexcept;
  // Creates scatter cross section from a scatter matrix by summing each column
  static OneDimensional
  CreateScatterXS(const pugi::xml_node& scatter_node);
  // Helper function for reaction cross section construction. Throws exception
  // on incorrect number of entries.
  static ReactionsMap
  CreateReactions(const pugi::xml_node& particle_node);
  // Helper function for creating total cross section by summing all reactions
  static OneDimensional CreateTotalXS(
      const pugi::xml_node& particle_node,
      const ReactionsMap& reactions) noexcept;
  // Captures the particle, killing it
  void Capture(State& s) const noexcept;
  // Scatters the particle and updates its Group and Direction
  void Scatter(State& s) const noexcept;
  // Fissions the Nuclide and produces secondaries
  void Fission(State& s) const noexcept;
  // Average number of secondary particles produced per fission
  const std::optional<OneDimensional> nubar;
  // Outgoing energy distribution of fission neutrons
  const std::optional<NormalizedTwoDimensional> chi;
  // Probabilities for scattering from incoming energy to outgoing energy
  const std::optional<NormalizedTwoDimensional> scatter_probs;
  // Cross section data for each Reactions
  const ReactionsMap reactions;
  // Total cross section obtained by summing all reactions
  const OneDimensional total;
  // Number of groups
  const Group max_group;
};
