#pragma once

#include "BasicTypes.hpp"
#include "HDF5DataSet.hpp"
#include "Interaction.hpp"
#include "Reaction.hpp"
#include "pugixml.hpp"

#include <filesystem>
#include <map>
#include <optional>

class Particle;

/// @brief Contains cross sections which are indexed by continuous energy values
class Continuous : public Interaction {
public:
  /// @brief Constructs continuous energy nuclear data from a particle node of
  ///        an XML document
  Continuous(const pugi::xml_node& particle_node);
  /// @brief Returns the total cross section for a given Particle
  MicroscopicCrossSection GetTotal(const Particle& p) const noexcept override;
  /// @brief Returns the cross section for a given Particle and Reaction
  MicroscopicCrossSection
  GetReaction(const Particle& p, const Reaction r) const noexcept override;
  /// @brief Returns the average fission neutron yield for a given Particle
  Real GetNuBar(const Particle& p) const noexcept override;
  /// @brief Interact with a Particle, updating its state
  void Interact(Particle& p) const noexcept override;

private:
  // Continuously maps elements from a domain to a range, provided a limited
  // set of points, by interpolating
  template <typename Key, typename T> class Map {
  public:
    // Type used to store elements internally
    using elements_type = std::map<Key, T>;
    // Constructs Continuous::Map by assigning elements directly
    Map(elements_type&& other);
    // Returns a const reference to the value at a given key
    const T& at(const Key k) const noexcept;

  protected:
    // This class essentially wraps an STL container
    const elements_type elements;
  };
  // Like Map, but stores elements as the CDF of some random variable. Stores
  // CDF values as keys so that std::map::upper_bound() can be used.
  template <typename T> class CDF : public Map<Real, T> {
  public:
    // Constructs a CDF from a std::map
    CDF(typename Map<Real, T>::elements_type&& other);
    // Samples a value from the CDF and returns the sampled key
    const T& Sample(RNG& rng) const noexcept;
  };
  // Continuous energy cross sections
  using CE_XS = Map<ContinuousEnergy, MicroscopicCrossSection>;
  // Helper function for constructing CE_XS from JANIS Web data file
  static CE_XS::elements_type
  ReadJanisWeb(const std::filesystem::path& datapath);
  // Helper function for constructing CDF<ContinuousEnergy> from JANIS Web data
  // file
  static CDF<ContinuousEnergy>::elements_type
  ReadJanisWebCDF(const std::filesystem::path& datapath);
  // Helper function for constructing thermal scattering data S(a,b,T) if found
  static std::optional<HDF5DataSet>
  ReadPandasSAB(const pugi::xml_node& tsl_node);
  // Helper function for reaction cross section construction
  static std::map<Reaction, CE_XS>
  CreateReactions(const pugi::xml_node& particle_node);
  // Captures the Particle, killing it
  void Capture(Particle& p) const noexcept;
  // Scatters the Particle and updates its ContinuousEnergy and Direction
  void Scatter(Particle& p) const noexcept;
  /// @brief Fissions the Nuclide and produces secondaries
  void Fission(Particle& p) const noexcept;
  // Average number of secondary particles produced per fission
  const std::optional<Map<ContinuousEnergy, Real>> nubar;
  // Outgoing energy distribution of fission neutrons
  const std::optional<CDF<ContinuousEnergy>> chi;
  // Neutron thermal scattering law S(a,b,T)
  const std::optional<HDF5DataSet> sab;
  // Cross section data for each Reaction
  const std::map<Reaction, CE_XS> reactions;
  // Total cross section provided in nuclear data files
  const CE_XS total;
};
