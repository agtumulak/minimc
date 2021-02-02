#include "Transport.hpp"

#include "Constants.hpp"

#include <algorithm>
#include <iterator>

TransportOutcome
TransportAndBank(const TransportOutcome& input, const World& w) {
  Estimator e = input.estimator;
  Bank b;
  for (auto p : input.banked) {
    RNG rng{p.seed};
    while (p.IsAlive()) {
      const auto collision = p.GetCell().SampleCollisionDistance(rng, p);
      const auto [nearest_surface, surface_crossing] =
          p.GetCell().NearestSurface(p.GetPosition(), p.GetDirection());
      if (collision < surface_crossing) {
        e.at(Estimator::Event::collision)++;
        p.Stream(collision);
        switch (const auto& nuclide =
                    p.GetCell().material->SampleNuclide(rng, p);
                nuclide.SampleReaction(rng, p)) {
        case NuclearData::Reaction::capture:
          e.at(Estimator::Event::capture)++;
          p.Kill();
          break;
        case NuclearData::Reaction::scatter:
          e.at(Estimator::Event::scatter)++;
          nuclide.Scatter(rng, p);
          break;
        case NuclearData::Reaction::fission:
          e.at(Estimator::Event::fission)++;
          auto fission_yield{nuclide.Fission(rng, p)};
          std::move(
              fission_yield.begin(), fission_yield.end(),
              std::back_insert_iterator(b));
          break;
        }
      }
      else {
        e.at(Estimator::Event::surface_crossing)++;
        p.Stream(surface_crossing + constants::nudge);
        p.SetCell(w.FindCellContaining(p.GetPosition()));
        if (!p.GetCell().material) {
          p.Kill();
        }
      }
    }
  }
  return TransportOutcome{e, b};
}
