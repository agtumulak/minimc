#pragma once

#include <memory>

namespace pugi {
class xml_node;
}
class History;
class World;

/// @brief Interface for a method which generates a full realization of a
///        History from an initial State
class TransportMethod {
public:
  /// @brief Factory method to create new TransportMethod from an XML document
  ///        and World
  /// @details This class is meant to compose the various tracking techniques
  ///          (surface tracking, delta tracking, etc...) and variance
  ///          reduction techniques (splitting, rouletting, forced collision)
  ///          by assigning such tasks to delegates. World is passed to
  ///          determine if the requested TransportMethod is @ref
  ///          transport_methods "supported".
  /// @returns A `std::unique_ptr` to the constructed TransportMethod (C++ Core
  ///          Guidelines R.30)
  static std::unique_ptr<const TransportMethod>
  Create(const pugi::xml_node& root, const World& world);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~TransportMethod() noexcept;
  /// @brief Samples each subsequent State in the History
  /// @param h History which is to be completed
  /// @param w Transition probabilities between each State are parameterized by
  ///          the geometric and material properties of a World
  virtual void Transport(History& h, const World& w) const noexcept = 0;
};

/// @brief Loops over each CSGSurface in the current Cell to find the next
///        surface crossing
class SurfaceTracking : public TransportMethod {
public:
  /// @brief Implements surface tracking
  virtual void Transport(History& h, const World& w) const noexcept override;
};

/// @brief Tracks a particle using delta tracking within a Cell. Surface
///        tracking is used across different Cells.
class CellDeltaTracking : public TransportMethod {
public:
  /// @brief Implements cell delta tracking
  virtual void Transport(History& h, const World& w) const noexcept override;
};
