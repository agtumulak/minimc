#include "Source.hpp"

#include "World.hpp"

Source::Source() noexcept {}

Particle Source::Sample(std::minstd_rand& rng, const World& w) const noexcept {
  Particle p{};
  p.SetDirectionIsotropic(rng);
  p.SetCell(w.FindCellContaining(p.GetPosition()));
  return p;
}
