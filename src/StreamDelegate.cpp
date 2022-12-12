#include "StreamDelegate.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "Estimator/Proxy.hpp"
#include "Estimator/Visitor.hpp"
#include "Material.hpp"
#include "Particle.hpp"
#include "Perturbation/IndirectEffect/IndirectEffect.hpp"
#include "Perturbation/IndirectEffect/Visitor.hpp"
#include "Point.hpp"
#include "Reaction.hpp"
#include "World.hpp"
#include "pugixml.hpp"

#include <cassert>
#include <map>
#include <random>
#include <stdexcept>
#include <string>

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
    Particle& p, std::vector<Estimator::Proxy>& estimator_proxies,
    const Real distance) const noexcept {
  // move the Particle to new position
  p.SetPosition(p.GetPosition() + p.GetDirection() * distance);
  // update indirect effects
  for (auto& indirect_effect : p.indirect_effects) {
    indirect_effect->Visit(
        *GetStreamWithinCellIndirectEffectVisitor(p, distance));
  }
  // score Estimator objects
  for (auto& estimator_proxy : estimator_proxies) {
    estimator_proxy.Visit(*GetStreamWithinCellEstimatorVisitor(p, distance));
  }
}

void StreamDelegate::CrossSurface(
    Particle& p, std::vector<Estimator::Proxy>& estimator_proxies,
    const World& w, const CSGSurface& surface) const noexcept {
  const auto& new_cell = w.FindCellContaining(p.GetPosition());
  p.SetCell(new_cell);
  // update Estimator objects
  for (auto& estimator_proxy : estimator_proxies) {
    estimator_proxy.Visit(*GetCrossSurfaceEstimatorVisitor(p, surface));
  }
  // if Particle leaked, update Reaction and indirect effects
  if (new_cell.IsVoid()) {
    p.reaction = Reaction::leak;
  }
}

// SurfaceTracking

//// public

void SurfaceTracking::Stream(
    Particle& p, std::vector<Estimator::Proxy>& estimator_proxies,
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

std::unique_ptr<const Perturbation::IndirectEffect::Visitor>
SurfaceTracking::GetStreamWithinCellIndirectEffectVisitor(
    const Particle& p, const Real distance) const noexcept {
  // construct visitor for getting perturbation indirect effect
  class Visitor : public Perturbation::IndirectEffect::Visitor {
  public:
    Visitor(const Particle& p, Real distance) : p{p}, distance{distance} {}
    void Visit(Perturbation::IndirectEffect::TotalCrossSection&
                   total_xs_perturbation) const noexcept final {
      const auto& afracs = p.GetCell().material->afracs;
      if (afracs.find(total_xs_perturbation.nuclide) != afracs.cend()) {
        total_xs_perturbation.indirect_effects.front() +=
            1 / p.GetCell().material->GetMicroscopicTotal(p) - distance;
      }
    }

  private:
    const Particle& p;
    const Real distance;
  };
  return std::make_unique<const Visitor>(p, distance);
}

std::unique_ptr<const Estimator::Visitor>
SurfaceTracking::GetStreamWithinCellEstimatorVisitor(
    const Particle& p, const Real) const noexcept {
  class Visitor : public Estimator::Visitor {
  public:
    Visitor(const Particle& p) noexcept : Estimator::Visitor{p} {};
    Score Visit(const Estimator::Current&) const noexcept final { return 0; }
  };
  return std::make_unique<const Visitor>(p);
}

std::unique_ptr<const Estimator::Visitor>
SurfaceTracking::GetCrossSurfaceEstimatorVisitor(
    const Particle& p, const CSGSurface& s) const noexcept {
  class Visitor : public Estimator::Visitor {
  public:
    Visitor(const Particle& p, const CSGSurface& s)
        : Estimator::Visitor{p}, s{s} {}
    Score
    Visit(const Estimator::Current& current_estimator) const noexcept final {
      return current_estimator.surface.get() == &s ? 1 : 0;
    }

  private:
    const CSGSurface& s;
  };
  return std::make_unique<const Visitor>(p, s);
}

// CellDetltaTracking

//// public

void CellDeltaTracking::Stream(
    Particle& p, std::vector<Estimator::Proxy>& estimator_proxies,
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

std::unique_ptr<const Perturbation::IndirectEffect::Visitor>
CellDeltaTracking::GetStreamWithinCellIndirectEffectVisitor(
    const Particle& p, const Real) const noexcept {
  // construct visitor for getting perturbation indirect effect
  class Visitor : public Perturbation::IndirectEffect::Visitor {
  public:
    Visitor(const Particle& p) : p{p} {}
    void Visit(Perturbation::IndirectEffect::TotalCrossSection&
                   total_xs_perturbation) const noexcept final {
      const auto& afracs = p.GetCell().material->afracs;
      if (afracs.find(total_xs_perturbation.nuclide) != afracs.cend()) {
        total_xs_perturbation.indirect_effects.front() +=
            1 / p.GetCell().material->GetMicroscopicTotal(p);
      }
    }

  private:
    const Particle& p;
  };
  return std::make_unique<const Visitor>(p);
}

std::unique_ptr<const Estimator::Visitor>
CellDeltaTracking::GetStreamWithinCellEstimatorVisitor(
    const Particle& p, const Real) const noexcept {
  class Visitor : public Estimator::Visitor {
  public:
    Visitor(const Particle& p) noexcept : Estimator::Visitor{p} {};
    Score Visit(const Estimator::Current&) const noexcept final { return 0; }
  };
  return std::make_unique<const Visitor>(p);
}

std::unique_ptr<const Estimator::Visitor>
CellDeltaTracking::GetCrossSurfaceEstimatorVisitor(
    const Particle& p, const CSGSurface& s) const noexcept {
  class Visitor : public Estimator::Visitor {
  public:
    Visitor(const Particle& p, const CSGSurface& s)
        : Estimator::Visitor{p}, s{s} {}
    Score
    Visit(const Estimator::Current& current_estimator) const noexcept final {
      return current_estimator.surface.get() == &s ? 1 : 0;
    }

  private:
    const CSGSurface& s;
  };
  return std::make_unique<const Visitor>(p, s);
}
