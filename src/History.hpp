#pragma once

#include "BasicTypes.hpp"
#include "State.hpp"

#include <memory>
#include <vector>

class Cell;
class CSGSurface;
class Source;
class World;

/// @brief Contains all information required to uniquely specify a history and
///        efficiently produce an Estimator score
class History {
public:
  /// @brief Constructs a History with an initial State sampled from a Source
  History(const World& w, const Source& source, RNG::result_type seed) noexcept;
  /// @brief Streams the particle across a CSGSurface into a new Cell
  void CrossSurface(
      std::shared_ptr<const CSGSurface> surface, const Real distance,
      const Cell& cell) noexcept;
  /// @brief Collides the particle at some distance along its direction of
  ///        flight without crossing a CSGSurface
  void CollideWithinCell(const Real distance) noexcept;
  /// @brief Streams the particle along its direction of flight for a given
  ///        distance by updating the current State
  /// @details Does not update most recent CSGSurface or current Cell
  void StreamWithinCell(const Real distance) noexcept;
  /// @brief Returns a const reference to the most recent State
  /// @details For a Markovian process, the most recent State is all that is
  ///          required to sample the next State
  const State& GetState() const noexcept;
  /// @brief Returns a reference to the RNG of the most recent State
  RNG& GetRNG() noexcept;
  /// @brief Returns a const iterator to the first State in the History
  std::vector<State>::const_iterator begin() const noexcept;
  /// @brief Returns a const iterator to the element following the last State
  ///        in the History
  std::vector<State>::const_iterator end() const noexcept;

private:
  // States previously occupied
  std::vector<State> history;
};
