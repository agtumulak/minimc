#pragma once

#include "BasicTypes.hpp"
#include "CSGSurface.hpp"

#include <memory>
#include <vector>

namespace Estimator {
class Proxy;
class Visitor;
} // namespace Estimator
namespace Perturbation {
namespace IndirectEffect {
class Visitor;
} // namespace IndirectEffect
} // namespace Perturbation
namespace pugi {
class xml_node;
}
class Particle;
class World;

/// @brief Models the streaming of a Particle through a World
/// @details Composition over inheritance interface for Stream() method
class StreamDelegate {
public:
  /// @brief Factory method to create new StreamDelegate from an XML document
  ///        and World
  /// @details World is passed to determine if the requested StreamDelegate
  ///          is @ref stream_delegates "supported".
  /// @returns A `std::unique_ptr` to the constructed StreamDelegate (C++ Core
  ///          Guidelines R.30)
  static std::unique_ptr<const StreamDelegate>
  Create(const pugi::xml_node& root, const World& world);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~StreamDelegate() noexcept;
  /// @brief Interface for scoring track length estimators and updating
  ///        indirect effects
  /// @details Transports a Particle up to (but not including) a collision
  virtual void Stream(
      Particle& p, std::vector<Estimator::Proxy>& estimator_proxies,
      const World& w) const noexcept = 0;

protected:
  /// @brief Interface for streaming a Particle within a Cell once a distance
  ///        has been sampled
  void StreamWithinCell(
      Particle& p, std::vector<Estimator::Proxy>& estimator_proxies,
      const Real distance) const noexcept;
  /// @brief Interface for moving a Particle across a CSGSurface while streaming
  /// @todo Change nudge to be in direction normal to surface crossed in order
  ///       to reduce minimum nudge size and/or avoid "grazing" intersections.
  void CrossSurface(
      Particle& p, std::vector<Estimator::Proxy>& estimator_proxies,
      const World& w, const CSGSurface& surface) const noexcept;

private:
  // Interface for updating each indirect effect after streaming within the
  // current Cell. Reference to const Particle ensures that order each
  // indirect effect is visited does not matter.
  virtual std::unique_ptr<const Perturbation::IndirectEffect::Visitor>
  GetStreamWithinCellIndirectEffectVisitor(
      const Particle& p, const Real distance) const noexcept = 0;
  // Interface for getting the score that streaming within a Cell would
  // contribute to an estimator
  virtual std::unique_ptr<const Estimator::Visitor>
  GetStreamWithinCellEstimatorVisitor(
      const Particle& p, const Real distance) const noexcept = 0;
  // Interface for getting the score that streaming across a CSGSurface would
  // contribute to an estimator
  virtual std::unique_ptr<const Estimator::Visitor>
  GetCrossSurfaceEstimatorVisitor(
      const Particle& p, const CSGSurface& s) const noexcept = 0;
};

/// @brief Loops over each CSGSurface in the current Cell to find the next
///        surface crossing
class SurfaceTracking : public StreamDelegate {
public:
  void Stream(
      Particle& p, std::vector<Estimator::Proxy>& estimator_proxies,
      const World& w) const noexcept;

private:
  std::unique_ptr<const Perturbation::IndirectEffect::Visitor>
  GetStreamWithinCellIndirectEffectVisitor(
      const Particle& p, const Real distance) const noexcept final;

  std::unique_ptr<const Estimator::Visitor> GetStreamWithinCellEstimatorVisitor(
      const Particle& p, const Real distance) const noexcept final;

  std::unique_ptr<const Estimator::Visitor> GetCrossSurfaceEstimatorVisitor(
      const Particle& p, const CSGSurface& s) const noexcept final;
};

/// @brief Tracks a particle using delta tracking within a Cell. CSGSurface
///        tracking is used across different Cells.
class CellDeltaTracking : public StreamDelegate {
public:
  void Stream(
      Particle& p, std::vector<Estimator::Proxy>& estimator_proxies,
      const World& w) const noexcept override;

private:
  std::unique_ptr<const Perturbation::IndirectEffect::Visitor>
  GetStreamWithinCellIndirectEffectVisitor(
      const Particle& p, const Real distance) const noexcept final;

  std::unique_ptr<const Estimator::Visitor> GetStreamWithinCellEstimatorVisitor(
      const Particle& p, const Real distance) const noexcept final;

  std::unique_ptr<const Estimator::Visitor> GetCrossSurfaceEstimatorVisitor(
      const Particle& p, const CSGSurface& s) const noexcept final;
};
