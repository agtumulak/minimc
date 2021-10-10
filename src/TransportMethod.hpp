#pragma once

#include "Bank.hpp"
#include "pugixml.hpp"

#include <memory>

class Particle;
class EstimatorSet;
class World;

/// @brief Performs the transport of a Particle after it is born up until its
///        death
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
  /// @brief Moves a Particle through its states until it dies.
  /// @param p Particle to transport
  /// @param e EstimatorSet to score upon, typically owned by a single thread
  /// @param w World the Particle transports within
  virtual Bank
  Transport(Particle& p, EstimatorSet& e, const World& w) const noexcept = 0;
};

/// @brief Loops over each CSGSurface in the current Cell to find the next
///        surface crossing
class SurfaceTracking : public TransportMethod {
public:
  /// @brief Implements delta tracking
  Bank Transport(
      Particle& p, EstimatorSet& e, const World& w) const noexcept override;
};

/// @brief Tracks a particle using delta tracking within a Cell. Surface
///        tracking is used across different Cells.
class CellDeltaTracking : public TransportMethod {
public:
  /// @brief Implements cell delta tracking
  Bank Transport(
      Particle& p, EstimatorSet& e, const World& w) const noexcept override;
};
