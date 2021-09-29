#include "TransportMethod.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "Estimator.hpp"
#include "Material.hpp"
#include "Nuclide.hpp"
#include "Reaction.hpp"
#include "World.hpp"

#include <cassert>
#include <random>
#include <stdexcept>
#include <string>

// TransportMethod::Outcome

//// public

TransportMethod::Outcome&
TransportMethod::Outcome::operator+=(Outcome&& rhs) noexcept {
  estimator += rhs.estimator;
  banked.splice(banked.begin(), rhs.banked);
  return *this;
}

// TransportMethod

//// public

std::unique_ptr<const TransportMethod>
TransportMethod::Create(const pugi::xml_node& root, const World& world) {
  std::unique_ptr<const TransportMethod>
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

TransportMethod::~TransportMethod() noexcept {}

// SurfaceTracking

// SurfaceTracking

TransportMethod::Outcome
SurfaceTracking::Transport(Particle& p, const World& w) const noexcept {
  Outcome result;
  p.SetCell(w.FindCellContaining(p.GetPosition()));
  while (p.alive) {
    const auto distance_to_collision = std::exponential_distribution{
        p.GetCell().material->number_density *
        p.GetCell().material->GetMicroscopicTotal(p)}(p.rng);
    const auto [nearest_surface, distance_to_surface_crossing] =
        p.GetCell().NearestSurface(p.GetPosition(), p.GetDirection());
    if (distance_to_collision < distance_to_surface_crossing) {
      result.estimator.at(Estimator::Event::collision) += 1;
      p.Stream(distance_to_collision);
      const auto& nuclide = p.SampleNuclide();
      result.estimator.at(Estimator::Event::implicit_fission) +=
          nuclide.GetNuBar(p) * nuclide.GetReaction(p, Reaction::fission) /
          nuclide.GetTotal(p);
      nuclide.Interact(p);
    }
    else {
      result.estimator.at(Estimator::Event::surface_crossing) += 1;
      p.Stream(distance_to_surface_crossing + constants::nudge);
      p.SetCell(w.FindCellContaining(p.GetPosition()));
      if (!p.GetCell().material) {
        p.Kill();
      }
    }
  }
  return result;
}

// CellDeltaTracking

// TODO: Remove embarrasing amount of code duplication between this and
// SurfaceTracking
TransportMethod::Outcome
CellDeltaTracking::Transport(Particle& p, const World& w) const noexcept {
  Outcome result;
  p.SetCell(w.FindCellContaining(p.GetPosition()));
  while (p.alive) {
    const auto distance_to_collision = std::exponential_distribution{
        p.GetCell().material->number_density *
        p.GetCell().material->GetMicroscopicMajorant(p)}(p.rng);
    const auto [nearest_surface, distance_to_surface_crossing] =
        p.GetCell().NearestSurface(p.GetPosition(), p.GetDirection());
    if (distance_to_surface_crossing < distance_to_collision) {
      // a surface crossing occurs
      p.Stream(distance_to_surface_crossing + constants::nudge);
      result.estimator.at(Estimator::Event::surface_crossing) += 1;
      p.SetCell(w.FindCellContaining(p.GetPosition()));
      if (!p.GetCell().material){
        p.Kill();
      }
    }
    else if (std::bernoulli_distribution{
                 p.GetCell().material->GetMicroscopicTotal(p) /
                 p.GetCell().material->GetMicroscopicMajorant(p)}(p.rng)) {
      // a real collision occurs
      p.Stream(distance_to_collision);
      result.estimator.at(Estimator::Event::collision) += 1;
      const auto& nuclide = p.SampleNuclide();
      result.estimator.at(Estimator::Event::implicit_fission) +=
          nuclide.GetNuBar(p) * nuclide.GetReaction(p, Reaction::fission) /
          nuclide.GetTotal(p);
      nuclide.Interact(p);
    }
    else {
      // a virtual collision occurs, do nothing
      p.Stream(distance_to_collision);
    }
  }
  return result;
}
