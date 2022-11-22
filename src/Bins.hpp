#pragma once

#include "BasicTypes.hpp"
#include "Point.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <iterator>
#include <memory>
#include <optional>
#include <vector>

namespace pugi {
class xml_node;
}
class Particle;

/// @brief Partitions @f$ \mathbb{R} @f$ according to user-specified boundaries
/// @details A Bins defines boundaries @f$ a_{1}, \ldots, a_{N} @f$ which
///          define @f$ N + 1 @f$ bins @f$ b_{0}, \ldots, b_{N} @f$. Bin @f$ i
///          @f$ covers the range @f$ [ a_{i}, a_{i+1} ) @f$. Bin @f$ b_{0} @f$
///          has an implied lower boundary of @f$ a_{0} = -\infty @f$ and
///          bin @f$ b_{N} @f$ has an implied upper boundary @f$ a_{N+1} =
///          +\infty @f$.
class Bins {
public:
  /// @brief Factory method to create a new Bins from a BinTypeBase node of an
  ///        XML document
  static std::unique_ptr<const Bins> Create(const pugi::xml_node& bins_node);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Bins() noexcept;
  /// @brief Abstract interface for getting then number of bins @f$ N + 1 @f$
  /// @details Must return the actual size of container required to store all
  ///          bins
  virtual size_t size() const noexcept = 0;
  /// @brief Abstract interface for getting the bin index of a given value
  virtual BinIndex GetIndex(const Real& v) const noexcept = 0;
  /// @brief Returns the contents of Bins as a string suitable for printing
  virtual std::string to_string() const noexcept = 0;
};

/// @brief A Bin where there are no bins
class NoBins : public Bins {
public:
  /// @brief A NoBins only has one element
  size_t size() const noexcept final;
  /// @brief A NoBins always returns the one and only bin index, 0
  BinIndex GetIndex(const Real&) const noexcept final;
  /// @brief Prints string "none" to indicate absence of bins
  std::string to_string() const noexcept final;
};

/// @brief A Bin with equally-spaced bin boundaries
class LinspaceBins : public Bins {
public:
  /// @brief Constructs a LinspaceBins from a `linspace` node of an XML
  ///        document
  /// @details The `bins` attribute is interpreted to mean the number of bins
  ///          of <em>finite</em> size. Two additional bins of infinte-width
  ///          are automatically included.
  LinspaceBins(const pugi::xml_node& linspace_node);
  /// @brief Size required to store all bins
  size_t size() const noexcept final;
  /// @brief The equally spaced bin structure can take advantage of arithmetic
  ///        that speeds up bin index lookup
  BinIndex GetIndex(const Real& v) const noexcept final;
  /// @brief Prints space-separated array of bin boundaries
  std::string to_string() const noexcept final;

private:
  // total number of bins
  const size_t n_bins;
  // lowest finite bin boundary
  const Real lower_bound;
  // highest finite bin boundary
  const Real upper_bound;
  // bin width
  const Real bin_width;
};

/// @brief A Bin with logarithmically-spaced bin boundaries
class LogspaceBins : public Bins {
public:
  /// @brief Constructs a LogspaceBins from a `logspace` node of an XML
  ///        document
  /// @details The `bins` attribute is interpreted to mean the number of bins
  ///          of <em>finite</em> size. Two additional bins of infinte-width
  ///          are automatically included.
  LogspaceBins(const pugi::xml_node& logspace_node);
  /// @brief Size required to store all bins
  size_t size() const noexcept final;
  /// @brief The logarithmically spaced bin structure can take advantage of
  ///        arithmetic that speeds up bin index lookup
  BinIndex GetIndex(const Real& v) const noexcept final;
  /// @brief Prints space-separated array of bin boundaries
  std::string to_string() const noexcept final;

private:
  // number of bins; equal to number of bin boundaries
  const size_t n_bins;
  // base of the log space
  const Real base;
  // lowest finite bin boundary in log space
  const Real log_lower_bound;
  // highest finite bin boundary in log space
  const Real log_upper_bound;
  // bin width in log space
  const Real log_bin_width;
};

/// @brief A Bin with explicitly specified bin boundaries
class BoundaryBins : public Bins {
public:
  /// @brief Constructs a BoundaryBins from a `boundaries` node of an XML
  ///        document
  BoundaryBins(const pugi::xml_node& boundaries_node);
  /// @brief Returns the number of elements
  size_t size() const noexcept final;
  /// @brief Return an index to the bin that would contain a given value
  BinIndex GetIndex(const Real& v) const noexcept final;
  /// @brief Prints space-separated array of bin boundaries
  std::string to_string() const noexcept final;

private:
  // bin boundaries
  const std::vector<Real> boundaries;
};

/// @brief Partitions all @ref estimators_phase_space "phase space" into
///        distinct bins.
class ParticleBins {
public:
  /// @brief Constructs a ParticleBins from a `bins` node of an XML document
  ParticleBins(const pugi::xml_node& bins_node) noexcept;
  /// @brief Total number of bins used
  size_t size() const noexcept;
  /// @brief Maps a point in @ref estimators_phase_space "phase space" @f$ X
  ///        @f$ to an index in @f$ \mathbb{Z}_{\geq 0} @f$
  BinIndex GetIndex(const Particle& p) const noexcept;
  /// @brief Returns a string suitable for printing
  std::string to_string() const noexcept;

private:
  // Helper visitor to convert Energy to Real
  struct VisitEnergy {
    Real operator()(const ContinuousEnergy& e) const noexcept { return e; }
    Real operator()(const Group& g) const noexcept { return g; }
  };
  // precomputes the stride for each dimension, base case
  template <typename T, typename U>
  static std::array<size_t, 2>
  ComputeStrides(const T&, const U& inner_bin) noexcept {
    return {inner_bin.size(), 1};
  }
  // precompute the strides for each dimension, general case
  template <typename T, typename U, typename... Args>
  static std::array<size_t, sizeof...(Args) + 2> ComputeStrides(
      const T&, const U& middle_bin, const Args&... inner_bins) noexcept {
    std::array<size_t, sizeof...(Args) + 2> strides;
    const auto inner_strides = ComputeStrides(middle_bin, inner_bins...);
    std::copy(
        inner_strides.cbegin(), inner_strides.cend(),
        std::next(strides.begin()));
    strides.front() = middle_bin.size() * inner_strides.front();
  }
  // reference direction against which direction cosines are computed
  const std::optional<Direction> direction;
  // cosine bins, uses direction member variable
  const std::unique_ptr<const Bins> cosine;
  // energy bins
  const std::unique_ptr<const Bins> energy;
  // contains the size of each element to facilitate constant-time lookup
  const std::array<size_t, 2> strides;
};
