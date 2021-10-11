#pragma once

#include "BasicTypes.hpp"
#include "pugixml.hpp"

#include <cstddef>
#include <iosfwd>
#include <memory>

/// @brief Partitions @f$ \mathbb{R} @f$ according to user-specified boundaries
/// @details A Bins defines boundaries @f$ a_{0}, \ldots, a_{N-1} @f$ which
///          define @f$ N @f$ bins. The @f$ i @f$-th bin contains all values
///          @f$ v @f$ which satisfy @f$ a_{i-1} \leq v < a_{i} @f$. The first
///          bin is implied to have a lower boundary of @f$ a_{-1} = -\infty
///          @f$. Values for which @f$ a_{N-1} \leq v @f$ throw an exception.
class Bins {
  /// @brief Writes contents of Bins to ostream; uses virtual friend function
  ///        idiom
  friend std::ostream& operator<<(std::ostream& os, const Bins& b) noexcept;
public:
  /// @brief Factory method to create a new Bins from a bins node of an XML
  ///        document
  static std::unique_ptr<const Bins>
  Create(const pugi::xml_node& bins_node) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Bins() noexcept;
  /// @brief Abstract interface for getting then number of bins
  virtual size_t size() const noexcept = 0;
  /// @brief Abstract interface for getting the bin index of a given value
  /// @exception std::out_of_range If the given value @f$ v @f$ and the last
  ///            bin boundary @f$ a_{N-1} @f$ satisfy @f$ a_{N-1} \leq v @f$,
  ///            then an exception is thrown.
  virtual size_t GetIndex(const Real& v) const = 0;

private:
  /// @brief Interface used by virtual friend function operator<<
  virtual void Print(std::ostream& os) const noexcept = 0;
};

/// @brief A Bin where there are no bins
class NoBins : public Bins {
public:
  /// @brief A NoBins only has one element
  size_t size() const noexcept override;
  /// @brief A NoBins always returns the one and only bin index, 0
  size_t GetIndex(const Real&) const noexcept override;

private:
  /// @brief Prints string "none" to indicate absence of bins
  void Print(std::ostream& os) const noexcept override;
};

/// @brief A Bin with equally-spaced bin boundaries
/// @todo Handle case where `min` > `max`
class LinspaceBins : public Bins {
public:
  /// @brief Constructs a LinspaceBins from a `linspace` node of an XML
  ///        document
  LinspaceBins(const pugi::xml_node& linspace_node) noexcept;
  /// @brief User-specified number of bins
  size_t size() const noexcept override;
  /// @brief The equally spaced bin structure can take advantage of arithmetic
  ///        that speeds up bin index lookup
  size_t GetIndex(const Real& v) const override;

private:
  void Print(std::ostream& os) const noexcept override;
  // number of bins; equal to number of bin boundaries
  const size_t n_bins;
  // lowest finite bin boundary
  const Real lower_bound;
  // bin width
  const Real bin_width;
};
