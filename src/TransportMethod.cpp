#include "TransportMethod.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "Estimator.hpp"
#include "Material.hpp"
#include "Nuclide.hpp"
#include "Particle.hpp"
#include "World.hpp"
#include "pugixml.hpp"

#include <cassert>
#include <iosfwd>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>

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

Bank SurfaceTracking::Transport(
    Particle& p, EstimatorSet& e, const World& w) const noexcept {
  Bank history_bank;
  p.SetCell(w.FindCellContaining(p.GetPosition()));
  while (p.IsAlive()) {
    const auto distance_to_collision = std::exponential_distribution{
        p.GetCell().material->number_density *
        p.GetCell().material->GetMicroscopicTotal(p)}(p.rng);
    const auto [nearest_surface, distance_to_surface_crossing] =
        p.GetCell().NearestSurface(p.GetPosition(), p.GetDirection());
    if (distance_to_collision < distance_to_surface_crossing) {
      p.Stream(distance_to_collision);
      const auto& nuclide = p.SampleNuclide();
      nuclide.Interact(p);
    }
    else {
      p.Stream(distance_to_surface_crossing + constants::nudge);
      p.SetCell(w.FindCellContaining(p.GetPosition()));
      p.current_surface = nearest_surface;
      p.event = p.GetCell().material ? Particle::Event::surface_cross
                                     : Particle::Event::leak;
    }
    e.Score(p);
  }
  return history_bank;
}

// CellDeltaTracking

// TODO: Remove embarrasing amount of code duplication between this and
// SurfaceTracking
Bank CellDeltaTracking::Transport(
    Particle& p, EstimatorSet& e, const World& w) const noexcept {
  Bank history_bank;
  p.SetCell(w.FindCellContaining(p.GetPosition()));
  while (p.IsAlive()) {
    const auto distance_to_collision = std::exponential_distribution{
        p.GetCell().material->number_density *
        p.GetCell().material->GetMicroscopicMajorant(p)}(p.rng);
    const auto [nearest_surface, distance_to_surface_crossing] =
        p.GetCell().NearestSurface(p.GetPosition(), p.GetDirection());
    if (distance_to_surface_crossing < distance_to_collision) {
      // a surface crossing occurs
      p.Stream(distance_to_surface_crossing + constants::nudge);
      p.SetCell(w.FindCellContaining(p.GetPosition()));
      p.current_surface = nearest_surface;
      p.event = p.GetCell().material ? Particle::Event::surface_cross
                                     : Particle::Event::leak;
    }
    else if (std::bernoulli_distribution{
                 p.GetCell().material->GetMicroscopicTotal(p) /
                 p.GetCell().material->GetMicroscopicMajorant(p)}(p.rng)) {
      // a real collision occurs
      p.Stream(distance_to_collision);
      const auto& nuclide = p.SampleNuclide();
      nuclide.Interact(p);
    }
    else {
      // a virtual collision occurs, do nothing
      p.Stream(distance_to_collision);
      p.event = Particle::Event::virtual_collision;
    }
    e.Score(p);
  }
  return history_bank;
}
