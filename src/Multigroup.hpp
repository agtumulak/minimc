#pragma once

#include "BasicTypes.hpp"
#include "NuclearData.hpp"
#include "Particle.hpp"
#include "pugixml.hpp"

#include <map>
#include <random>
#include <vector>

/// @brief Contains cross sections which are indexed by discrete energy groups
/// @details Groups are integers in `[1,G]`. Group `1` corresponds to the
///          highest energy. Group `G` corresponds to the lowest energy.
class Multigroup : public NuclearData {
public:
  /// @brief Constructs multigroup nuclear data from a particle node of an XML
  ///        document
  /// @details Creates cross sections for each Reaction. If no corresponding
  ///          reaction node is found, it is not added to the map.
  /// @param particle_node The requested particle node in the XML document
  /// @param G The largest group number expected. Corresponds to group with the
  ///        <em>lowest</em> energy.
  /// @exception std::runtime_error Number of entries is not consistent with
  ///            `G` parameter
  Multigroup(const pugi::xml_node& particle_node, const Group G);
  /// @brief Constructs multigroup nuclear data from an existing set of
  ///        Multigroup data and associated weights.
  Multigroup(const std::map<NuclearData, Real>& weights);
  /// @brief Returns the total cross section for a given Particle
  CrossSection GetTotal(const Particle& p) const noexcept override;
  /// @brief Scatters the Particle and updates its group and direction
  /// @exception std::runtime_error A suitable outgoing Group could not be
  ///            sampled.
  void Scatter(std::minstd_rand& rng, Particle& p) const override;
  /// @brief Samples a Reaction
  Reaction SampleReaction(
      std::minstd_rand& rng, const Particle& p) const noexcept override;

private:
  // Groupwise cross sections indexed by one Group
  class OneDimensional {
  private:
    // The STL container underlying a OneDimensional object
    using elements_type = std::vector<Real>;

  public:
    // Constructs Multigroup::OneDimensional from vector of Real
    OneDimensional(const elements_type& elements);
    // Returns a reference to the `g`-th group
    Real& at(const Group g);
    // Returns a const reference to the `g`-th group
    const Real& at(const Group g) const;
    // Returns an iterator corresponding to most energetic Group
    elements_type::const_iterator begin() const noexcept;
    // Returns an iterator corresponding to the Group following the least
    // energetic Group
    elements_type::const_iterator end() const noexcept;

  private:
    elements_type elements;
  };
  // Groupwise cross sections indexed by two Groups
  class TwoDimensional {
  private:
    // The STL container underlying a TwoDimensional object
    using elements_type = std::vector<OneDimensional>;

  public:
    // Constructs Multigroup::TwoDimensional from vector of vector of Real
    TwoDimensional(const elements_type& elements);
    // Returns a reference to the `g`-th group
    OneDimensional& at(const Group g);
    // Returns a const reference to the `g`-th group
    const OneDimensional& at(const Group g) const;
    // Returns an iterator corresponding to most energetic Group
    elements_type::const_iterator begin() const noexcept;
    // Returns an iterator corresponding to the Group following the least
    // energetic Group
    elements_type::const_iterator end() const noexcept;

  private:
    elements_type elements;
  };
  using ReactionsMap = std::map<NuclearData::Reaction, OneDimensional>;
  // Helper function to create scattering matrix. Throws exception on incorrect
  // number of entries.
  static TwoDimensional
  CreateScatterMatrix(const pugi::xml_node& particle_node, const Group G);
  // Helper function for reaction cross section construction. Throws exception
  // on incorrect number of entries.
  static ReactionsMap CreateReactions(
      const pugi::xml_node& particle_node, const Group G,
      const TwoDimensional& scatter_matirx);
  // Helper function for creating total cross section by summing all reactions
  static OneDimensional
  CreateTotalXS(const ReactionsMap& reactions, const Group G) noexcept;
  // Returns the Reaction cross section for a given Particle
  CrossSection GetReaction(const Particle& p, const Reaction r) const noexcept;
  // Returns the (unnormalized) outgoing group probabilities for a given
  // Particle
  const OneDimensional&
  GetOutgoingScatterProbs(const Particle& p) const noexcept;
  // Probabilities for scattering from incoming energy to outgoing energy
  const TwoDimensional scatter_matrix;
  // Cross section data for each Reactions
  const ReactionsMap reactions;
  // Total cross section obtained by summing all reactions
  const OneDimensional total;
};
