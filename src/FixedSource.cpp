#include "FixedSource.hpp"

#include "History.hpp"

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
  Estimator worker_estimator;
  while (const auto range = chunk_giver.Next()) {
    if (!range) {
      break;
    }
    for (auto h = range->first; h < range->second; h++) {
      auto source_particle{Sample(h)};
      worker_estimator += History{source_particle, world}.Transport().estimator;
    }
  }
  return worker_estimator;
}

//// private

Particle FixedSource::Sample(RNG::result_type history) const noexcept {
  RNG rng{history};
  auto p = source.Sample(rng);
  p.SetCell(world.FindCellContaining(p.GetPosition()));
  return p;
}
