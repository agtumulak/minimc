#include "TransportMethod.hpp"

#include "Cell.hpp"
#include "History.hpp"
#include "Material.hpp"
#include "State.hpp"
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

void SurfaceTracking::Transport(History& h, const World& w) const noexcept {
  while (h.GetState().IsAlive()) {
    const auto& s = h.GetState(); // for functions which accept const State
    const auto collision_distance = std::exponential_distribution{
        s.cell->material->number_density *
        s.cell->material->GetMicroscopicTotal(s)}(h.GetRNG());
    const auto [nearest_surface, surface_crossing_distance] =
        s.cell->NearestSurface(s.position, s.direction);
    if (surface_crossing_distance < collision_distance) {
      h.CrossSurface(
          nearest_surface, surface_crossing_distance,
          w.FindCellContaining(s.position));
    }
    else {
      h.CollideWithinCell(collision_distance);
    }
  }
}

// CellDeltaTracking

void CellDeltaTracking::Transport(History& h, const World& w) const noexcept {
  while (h.GetState().IsAlive()) {
    const auto& s = h.GetState();
    const auto collision_distance = std::exponential_distribution{
        s.cell->material->number_density *
        s.cell->material->GetMicroscopicMajorant(s)}(h.GetRNG());
    const auto [nearest_surface, surface_crossing_distance] =
        s.cell->NearestSurface(s.position, s.direction);
    if (surface_crossing_distance < collision_distance) {
      // a surface crossing occurs
      h.CrossSurface(
          nearest_surface, surface_crossing_distance,
          w.FindCellContaining(s.position));
    }
    else if (std::bernoulli_distribution{
                 s.cell->material->GetMicroscopicTotal(s) /
                 s.cell->material->GetMicroscopicMajorant(s)}(h.GetRNG())) {
      // a real collision occurs
      h.CollideWithinCell(collision_distance);
    }
    else {
      // a virtual collision occurs
      h.StreamWithinCell(collision_distance);
    }
  }
}
