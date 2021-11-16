#include "State.hpp"

#include "Cell.hpp"
#include "Material.hpp"
#include "Nuclide.hpp"
#include "Source.hpp"
#include "World.hpp"

#include <random>
#include <cassert>
#include <type_traits>

// State

//// public

State::State(
    const World& world, const Source& source, RNG::result_type seed) noexcept
    : rng{seed}, position{source.SamplePosition(rng)},
      direction{source.SampleDirection(rng)}, energy{source.SampleEnergy(rng)},
      cell{&(world.FindCellContaining(position))}, event{Event::birth},
      particle{source.SampleParticle(rng)} {}

State::State(
    const State& origin, std::shared_ptr<const CSGSurface> surface,
    const Real distance, const Cell& cell) noexcept
    : rng{origin.rng}, position{origin.position + origin.direction * distance},
      direction{origin.direction}, energy{origin.energy}, surface{surface},
      cell{&cell}, event{Event::surface_cross}, particle{origin.particle} {}

State::State(const State& origin, const Real distance) noexcept
    : rng{origin.rng}, position{origin.position + origin.direction * distance},
      direction{origin.direction}, energy{origin.energy}, cell{origin.cell},
      nuclide{SampleNuclide(*cell, rng)}, event{origin.event},
      particle{origin.particle} {
  // then interact with the sampled Nuclide
  nuclide->Interact(*this);
}

std::shared_ptr<const Nuclide>
State::SampleNuclide(const Cell& cell, RNG& rng) noexcept {
  const auto& material = *cell.material;
  const MicroscopicCrossSection threshold =
      material.GetMicroscopicTotal(*this) *
      std::uniform_real_distribution{}(rng);
  MicroscopicCrossSection accumulated = 0;
  for (const auto& [nuclide_ptr, afrac] : material.afracs) {
    accumulated += afrac * nuclide_ptr->GetTotal(*this);
    if (accumulated > threshold) {
      return nuclide_ptr;
    }
  }
  // This should never be reached since Material total cross section is
  // computed from constituent Nuclide total cross sections
  assert(false);
}

Real State::Sample() noexcept { return std::uniform_real_distribution{}(rng); }

bool State::IsAlive() const noexcept {
  return event != Event::capture && event != Event::leak &&
         event != Event::fission;
}
