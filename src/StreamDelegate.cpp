#include "StreamDelegate.hpp"

#include "Bins.hpp"
#include "Cell.hpp"
#include "Constants.hpp"
#include "Estimator.hpp"
#include "Material.hpp"
#include "Particle.hpp"
#include "Perturbation.hpp"
#include "Point.hpp"
#include "World.hpp"
#include "pugixml.hpp"

#include <cassert>
#include <map>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>

// StreamDelegate

//// public

std::unique_ptr<const StreamDelegate>
StreamDelegate::Create(const pugi::xml_node& root, const World& world) {
  std::unique_ptr<const StreamDelegate>
      transport_method{}; // hope this triggers NRVO
  const std::string tracking_type =
      root.child("general").child("tracking").child_value();
  // default to surface tracking
  if (tracking_type.empty() || tracking_type == "surface") {
    if (world.HasConstantTemperature()) {
      transport_method = std::make_unique<const SurfaceTracking>();
    }
    else {
      throw std::runtime_error(
          "Surface tracking with continuous global temperature not allowed");
    }
  }
  else if (tracking_type == "cell delta") {
    transport_method = std::make_unique<const CellDeltaTracking>();
  }
  else {
    assert(false); // this should have been caught by the validator
  }
  return transport_method;
}

StreamDelegate::~StreamDelegate() noexcept {}

void StreamDelegate::StreamWithinCell(
    Particle& p, std::vector<EstimatorProxy>& estimator_proxies,
    const Real distance) const noexcept {
  // move the Particle to new position
  p.SetPosition(p.GetPosition() + p.GetDirection() * distance);
  // update indirect effects
  for (auto& [perturbation_ptr, indirect_effect] : p.indirect_effects) {
    indirect_effect += perturbation_ptr->Visit(
        *GetStreamWithinCellIndirectEffectVisitor(p, distance));
  }
  // score Estimator objects
  for (auto& estimator_proxy : estimator_proxies) {
    estimator_proxy.Visit(*GetStreamWithinCellEstimatorVisitor(p, distance));
  }
}

void StreamDelegate::CrossSurface(
    Particle& p, std::vector<EstimatorProxy>& estimator_proxies, const World& w,
    const CSGSurface& surface) const noexcept {
  const auto& new_cell = w.FindCellContaining(p.GetPosition());
  p.SetCell(new_cell);
  // update Estimator objects
  for (auto& estimator_proxy : estimator_proxies) {
    estimator_proxy.Visit(*GetCrossSurfaceEstimatorVisitor(p, surface));
  }
  // if Particle leaked, update Particle::Event
  if (new_cell.IsVoid()){
    p.event = Particle::Event::leak;
  }
}

// SurfaceTracking

//// public

void SurfaceTracking::StreamToNextCollision(
    Particle& p, std::vector<EstimatorProxy>& estimator_proxies,
    const World& w) const noexcept {
  while (true) {
    const auto distance_to_collision = std::exponential_distribution{
        p.GetCell().material->number_density *
        p.GetCell().material->GetMicroscopicTotal(p)}(p.rng);
    const auto [nearest_surface, distance_to_surface_crossing] =
        p.GetCell().NearestSurface(p.GetPosition(), p.GetDirection());
    // check if collision within the Cell has occured
    if (distance_to_collision < distance_to_surface_crossing) {
      StreamWithinCell(p, estimator_proxies, distance_to_collision);
      return;
    }
    // collision did not occur so Particle has streamed to adjacent Cell
    StreamWithinCell(
        p, estimator_proxies, distance_to_surface_crossing + constants::nudge);
    CrossSurface(p, estimator_proxies, w, *nearest_surface);
    if (p.GetCell().IsVoid()) {
      return;
    }
  }
}

//// private

std::unique_ptr<const Perturbation::Visitor>
SurfaceTracking::GetStreamWithinCellIndirectEffectVisitor(
    const Particle& p, const Real distance) const noexcept {
  // construct visitor for getting Perturbation indirect effect
  class IndirectEffectVisitor : public Perturbation::Visitor {
  public:
    IndirectEffectVisitor(const Particle& p, Real distance)
        : p{p}, distance{distance} {}
    Real Visit(const TotalCrossSectionPerturbation& perturbation)
        const noexcept final {
      const auto& afracs = p.GetCell().material->afracs;
      if (afracs.find(perturbation.nuclide) == afracs.cend()) {
        // indirect effect is zero if perturbed Nuclide is not in Cell
        return 0;
      }
      else {
        return 1 / p.GetCell().material->GetMicroscopicTotal(p) - distance;
      }
      return 0;
    }

  private:
    const Particle& p;
    const Real distance;
  };
  return std::make_unique<const IndirectEffectVisitor>(p, distance);
}

std::unique_ptr<const Estimator::Visitor>
SurfaceTracking::GetStreamWithinCellEstimatorVisitor(
    const Particle&, const Real) const noexcept {
  class StreamWithinCellEstimatorVisitor : public Estimator::Visitor {
  public:
    Visitor::T Visit(const CurrentEstimator&) const noexcept final {
      return std::nullopt;
    }
  };
  return std::make_unique<const StreamWithinCellEstimatorVisitor>();
}

std::unique_ptr<const Estimator::Visitor>
SurfaceTracking::GetCrossSurfaceEstimatorVisitor(
    const Particle& p, const CSGSurface& s) const noexcept {
  class CrossSurfaceEstimatorVisitor : public Estimator::Visitor {
  public:
    CrossSurfaceEstimatorVisitor(const Particle& p, const CSGSurface& s)
        : p{p}, s{s} {}
    Visitor::T
    Visit(const CurrentEstimator& current_estimator) const noexcept final {
      return current_estimator.surface.get() == &s
                 ? Visitor::T{std::make_tuple(
                       current_estimator.bins->GetIndex(p), 1)}
                 : std::nullopt;
    }

  private:
    const Particle& p;
    const CSGSurface& s;
  };
  return std::make_unique<const CrossSurfaceEstimatorVisitor>(p, s);
}

// CellDetltaTracking

//// public

void CellDeltaTracking::StreamToNextCollision(
    Particle& p, std::vector<EstimatorProxy>& estimator_proxies,
    const World& w) const noexcept {
  while (true) {
    const auto distance_to_collision = std::exponential_distribution{
        p.GetCell().material->number_density *
        p.GetCell().material->GetMicroscopicMajorant(p)}(p.rng);
    const auto [nearest_surface, distance_to_surface_crossing] =
        p.GetCell().NearestSurface(p.GetPosition(), p.GetDirection());
    // check if collision within the Cell has occured
    if (distance_to_collision < distance_to_surface_crossing) {
      StreamWithinCell(p, estimator_proxies, distance_to_collision);
      // check if collision was a true collision
      if (std::bernoulli_distribution{
              p.GetCell().material->GetMicroscopicTotal(p) /
              p.GetCell().material->GetMicroscopicMajorant(p)}(p.rng)) {
        return;
      }
      // collision was virtual so we do nothing and sample next distance
      continue;
    }
    // collision did not occur so Particle has streamed to adjacent Cell
    StreamWithinCell(
        p, estimator_proxies, distance_to_surface_crossing + constants::nudge);
    CrossSurface(p, estimator_proxies, w, *nearest_surface);
    if (p.GetCell().IsVoid()) {
      return;
    }
  }
}

//// private

std::unique_ptr<const Perturbation::Visitor>
CellDeltaTracking::GetStreamWithinCellIndirectEffectVisitor(
    const Particle& p, const Real) const noexcept {
  // construct visitor for getting Perturbation indirect effect
  class IndirectEffectVisitor : public Perturbation::Visitor {
  public:
    IndirectEffectVisitor(const Particle& p) : p{p} {}
    Real Visit(const TotalCrossSectionPerturbation& perturbation)
        const noexcept final {
      const auto& afracs = p.GetCell().material->afracs;
      if (afracs.find(perturbation.nuclide) == afracs.cend()) {
        // indirect effect is zero if perturbed Nuclide is not in Cell
        return 0;
      }
      else {
        return 1 / p.GetCell().material->GetMicroscopicTotal(p);
      }
      return 0;
    }

  private:
    const Particle& p;
  };
  return std::make_unique<const IndirectEffectVisitor>(p);
}

std::unique_ptr<const Estimator::Visitor>
CellDeltaTracking::GetStreamWithinCellEstimatorVisitor(
    const Particle&, const Real) const noexcept {
  class StreamWithinCellEstimatorVisitor : public Estimator::Visitor {
  public:
    Visitor::T Visit(const CurrentEstimator&) const noexcept final {
      return std::nullopt;
    }
  };
  return std::make_unique<const StreamWithinCellEstimatorVisitor>();
}

std::unique_ptr<const Estimator::Visitor>
CellDeltaTracking::GetCrossSurfaceEstimatorVisitor(
    const Particle& p, const CSGSurface& s) const noexcept {
  class CrossSurfaceEstimatorVisitor : public Estimator::Visitor {
  public:
    CrossSurfaceEstimatorVisitor(const Particle& p, const CSGSurface& s)
        : p{p}, s{s} {}
    Visitor::T
    Visit(const CurrentEstimator& current_estimator) const noexcept final {
      return current_estimator.surface.get() == &s
                 ? Visitor::T{std::make_tuple(
                       current_estimator.bins->GetIndex(p), 1)}
                 : std::nullopt;
    }

  private:
    const Particle& p;
    const CSGSurface& s;
  };
  return std::make_unique<const CrossSurfaceEstimatorVisitor>(p, s);
}
