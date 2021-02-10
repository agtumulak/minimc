#include "Transport.hpp"

#include "Constants.hpp"

#include <algorithm>
#include <iterator>

// TransportOutcome

//// public

TransportOutcome&
TransportOutcome::operator+=(const TransportOutcome& rhs) noexcept {
  estimator += rhs.estimator;
  std::move(
      rhs.banked.cbegin(), rhs.banked.cend(),
      std::back_insert_iterator(banked));
  return *this;
}

// TransportAndBank

//// public

TransportOutcome TransportAndBank(Bank& source, const World& w) {
  Estimator e;
  Bank secondaries;
  // temp
  for (auto& p : source) {
    RNG rng{p.seed};
    p.SetCell(w.FindCellContaining(p.GetPosition()));
    while (p.IsAlive()) {
      const auto collision = p.GetCell().SampleCollisionDistance(rng, p);
      const auto [nearest_surface, surface_crossing] =
          p.GetCell().NearestSurface(p.GetPosition(), p.GetDirection());
      if (collision < surface_crossing) {
        e.at(Estimator::Event::collision) += 1;
        p.Stream(collision);
        const auto& nuclide = p.GetCell().material->SampleNuclide(rng, p);
        e.at(Estimator::Event::implicit_fission) +=
            nuclide.GetNuBar(p) * nuclide.GetFission(p) / nuclide.GetTotal(p);
        switch (nuclide.SampleReaction(rng, p)) {
        case NuclearData::Reaction::capture:
          e.at(Estimator::Event::capture) += 1;
          p.Kill();
          break;
        case NuclearData::Reaction::scatter:
          e.at(Estimator::Event::scatter) += 1;
          nuclide.Scatter(rng, p);
          break;
        case NuclearData::Reaction::fission:
          e.at(Estimator::Event::fission) += 1;
          auto fission_yield{nuclide.Fission(rng, p)};
          std::move(
              fission_yield.begin(), fission_yield.end(),
              std::back_insert_iterator(secondaries));
          break;
        }
      }
      else {
        e.at(Estimator::Event::surface_crossing) += 1;
        p.Stream(surface_crossing + constants::nudge);
        p.SetCell(w.FindCellContaining(p.GetPosition()));
        if (!p.GetCell().material) {
          p.Kill();
        }
      }
    }
  }
  return TransportOutcome{e, secondaries};
}
