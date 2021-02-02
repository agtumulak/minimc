#include "FixedSource.hpp"

#include "Constants.hpp"
#include "Transport.hpp"

#include <future>
#include <iostream>
#include <numeric>
#include <optional>

// FixedSource

//// public

FixedSource::FixedSource(const pugi::xml_node& root)
    : Driver{root}, source{root.child("problemtype").child("fixedsource")} {};

void FixedSource::Solve() {
  std::vector<std::future<Estimator>> results;
  std::cout << "Spawning " << std::to_string(threads) << " threads working on "
            << std::to_string(batchsize) << " histories split into "
            << std::to_string(
                   batchsize / chunksize + (batchsize % chunksize != 0))
            << " chunks... " << std::endl;
  for (size_t i = 0; i < threads; i++) {
    results.push_back(std::async(&FixedSource::StartWorker, this));
  }
  estimators = std::reduce(
      results.begin(), results.end(), Estimator{},
      [](auto& accumulated, auto& future) {
        return accumulated += future.get();
      });
  std::cout << estimators;
}

Estimator FixedSource::StartWorker() {
  TransportOutcome worker_outcome;
  while (const auto range = chunk_giver.Next()) {
    if (!range) {
      break;
    }
    for (auto h = range->first; h < range->second; h++) {
      worker_outcome.banked.push_back(Sample(h));
      while (!(worker_outcome.banked.empty())) {
        worker_outcome = TransportAndBank(worker_outcome, world);
      }
    }
  }
  return worker_outcome.estimator;
}

//// private

Particle FixedSource::Sample(RNG::result_type history) const noexcept {
  // avoid zero seed with +1
  RNG rng{(history + 1) * constants::seed_stride};
  auto p = source.Sample(rng);
  p.SetCell(world.FindCellContaining(p.GetPosition()));
  return p;
}
