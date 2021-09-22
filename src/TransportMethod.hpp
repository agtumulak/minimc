#pragma once

#include "Estimator.hpp"
#include "Particle.hpp"
#include "pugixml.hpp"

#include <list>
#include <memory>

class World;

/// @brief Performs the transport of a Particle after it is born up until its
///        death
class TransportMethod {
public:
  /// @brief The result of a Transport call
  struct Outcome {
    /// @brief Adds the result of another transport result to this result
    Outcome& operator+=(Outcome&& rhs) noexcept;
    /// @brief Estimators scored during transport
    Estimator estimator;
    /// @brief Secondary particles banked during transport
    std::list<Particle> banked;
  };
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
  virtual Outcome Transport(Particle& p, const World& w) const noexcept = 0;

};

/// @brief Loops over each CSGSurface in the current Cell to find the next
///        surface crossing
class SurfaceTracking : public TransportMethod {
public:
  /// @brief Implements delta tracking
  Outcome Transport(Particle& p, const World& w) const noexcept override;
};

/// @brief Tracks a particle using delta tracking within a Cell. Surface
///        tracking is used across different Cells.
class CellDeltaTracking : public TransportMethod {
public:
  /// @brief Implements cell delta tracking
  Outcome Transport(Particle& p, const World& w) const noexcept override;
};
