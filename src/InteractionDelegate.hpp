#pragma once

#include "BasicTypes.hpp"
#include "ContinuousEvaluation.hpp"
#include "ContinuousReaction.hpp"
#include "ThermalScattering.hpp"

#include <map>
#include <memory>
#include <optional>
#include <vector>

namespace Estimator {
class Proxy;
}
namespace pugi {
class xml_node;
}
class Nuclide;
class Particle;
enum class Reaction;

/// @brief Models the interaction between a Particle and a Nuclide
/// @details Composition over inheritance interface for Interact() method
class InteractionDelegate {
public:
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~InteractionDelegate() noexcept;
  /// @brief Returns reference to thermal neutron scattering law data
  /// @exception std::runtime_error Thermal neutron scattering law data
  ///                               undefined for multigroup physics
  virtual const std::optional<ThermalScattering>& GetTNSL() const = 0;
  /// @brief Returns the majorant cross section for a given Particle
  virtual MicroscopicCrossSection
  GetCellMajorant(const Particle& p) const noexcept = 0;
  /// @brief Returns the total cross section for a given Particle
  virtual MicroscopicCrossSection
  GetTotal(const Particle& p) const noexcept = 0;
  /// @brief Interact with a Particle, updating its state
  virtual void Interact(
      Particle& p,
      std::vector<Estimator::Proxy>& estimator_proxies) const noexcept = 0;
};

/// @brief Contains cross sections which are indexed by continuous energy values
class Continuous : public InteractionDelegate {
public:
  /// @brief Constructs continuous energy nuclear data from a particle node of
  ///        an XML document and a target Nuclide
  Continuous(const pugi::xml_node& particle_node, const Nuclide& target);
  /// @brief Returns reference to optional thermal neutron scattering law data
  const std::optional<ThermalScattering>& GetTNSL() const final;
  /// @brief Returns the largest cross section that may be found within the
  ///        current Cell
  /// @details Currently this is the cross section at the majorant temperature
  ///          in the Cell
  MicroscopicCrossSection
  GetCellMajorant(const Particle& p) const noexcept override;
  /// @brief Returns the total cross section for a given Particle
  /// @details This is not guaranteed to be consistent with the sum of all
  ///          mutually exclusive reactions. The total cross section is meant
  ///          to be a user-provided quantity to speed up calculations.
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept override;
  /// @brief Interact with a Particle, updating its state
  void Interact(Particle& p, std::vector<Estimator::Proxy>& estimator_proxies)
      const noexcept override;

private:
  // Returns true if any reaction modifies the total cross section even if the
  // total cross section was evaluated at the target temperature
  bool ReactionsModifyTotal(const Particle& p) const noexcept;
  // Thermal neutron scattering law S(a,b,T). Must be initialized before
  // `reactions` since ContinuousScatter uses a reference to `tnsl`
  const std::optional<ThermalScattering> tnsl;
  // Cross section data for mutually exclusive reactions
  const std::vector<std::unique_ptr<const ContinuousReaction>> reactions;
  // Total cross section provided in nuclear data files
  const ContinuousEvaluation total;
};

/// @brief Contains cross sections which are indexed by discrete energy groups
/// @details Groups are integers in `[1,G]`. Group `1` corresponds to the
///          highest energy. Group `G` corresponds to the lowest energy.
/// @todo Refactor Multigroup::reactions into vector instead of map as is now
///       done with Continuous::reactions
class Multigroup : public InteractionDelegate {
public:
  /// @brief Constructs multigroup nuclear data from a particle node of an XML
  ///        document
  /// @details Creates cross sections for each Reaction. If no corresponding
  ///          reaction node is found, it is not added to the map.
  /// @param particle_node The requested particle node in the XML document
  /// @exception std::runtime_error Number of entries is not consistent with
  ///            number of groups
  Multigroup(const pugi::xml_node& particle_node);
  /// @brief Throws an exception for multigroup physics
  const std::optional<ThermalScattering>& GetTNSL() const final;
  /// @brief Returns the total cross section for a given Particle, currently
  ///        the same as GetTotal().
  MicroscopicCrossSection
  GetCellMajorant(const Particle& p) const noexcept override;
  /// @brief Returns the total cross section for a given Particle
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept override;
  /// @brief Interact with a Particle, updating its state
  void Interact(Particle& p, std::vector<Estimator::Proxy>& estimator_proxies)
      const noexcept override;

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
  static OneDimensional CreateScatterXS(const pugi::xml_node& scatter_node);
  // Helper function for reaction cross section construction. Throws exception
  // on incorrect number of entries.
  static ReactionsMap CreateReactions(const pugi::xml_node& particle_node);
  // Helper function for creating total cross section by summing all reactions
  static OneDimensional CreateTotalXS(
      const pugi::xml_node& particle_node,
      const ReactionsMap& reactions) noexcept;
  // Captures the Particle, killing it
  void Capture(Particle& p, std::vector<Estimator::Proxy>& estimator_proxies)
      const noexcept;
  // Scatters the Particle and updates its Group and Direction
  void Scatter(Particle& p, std::vector<Estimator::Proxy>& estimator_proxies)
      const noexcept;
  // Fissions the Nuclide and produces secondaries
  void Fission(Particle& p, std::vector<Estimator::Proxy>& estimator_proxies)
      const noexcept;
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
