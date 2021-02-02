#include "History.hpp"

#include "Constants.hpp"

#include <algorithm>
#include <iterator>

// History

//// public

History::History(Particle& source, const World& world) : p{source}, w{world} {}

History::Outcome History::Transport() {
  Outcome result;
  while (p.IsAlive()) {
    const auto collision = p.GetCell().SampleCollisionDistance(rng, p);
    const auto [nearest_surface, surface_crossing] =
        p.GetCell().NearestSurface(p.GetPosition(), p.GetDirection());
    if (collision < surface_crossing) {
      result.estimator.at(Estimator::Event::collision)++;
      p.Stream(collision);
      switch (const auto& nuclide = p.GetCell().material->SampleNuclide(rng, p);
              nuclide.SampleReaction(rng, p)) {
      case NuclearData::Reaction::capture:
        result.estimator.at(Estimator::Event::capture)++;
        p.Kill();
        break;
      case NuclearData::Reaction::scatter:
        result.estimator.at(Estimator::Event::scatter)++;
        nuclide.Scatter(rng, p);
        break;
      case NuclearData::Reaction::fission:
        result.estimator.at(Estimator::Event::fission)++;
        auto fission_yield{nuclide.Fission(rng, p)};
        std::move(
            fission_yield.begin(), fission_yield.end(),
            std::back_insert_iterator(result.banked));
        break;
      }
    }
    else {
      result.estimator.at(Estimator::Event::surface_crossing)++;
      p.Stream(surface_crossing + constants::nudge);
      p.SetCell(w.FindCellContaining(p.GetPosition()));
      if (!p.GetCell().material) {
        p.Kill();
      }
    }
  }
  return result;
}
