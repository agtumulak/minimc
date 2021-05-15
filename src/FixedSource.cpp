#include "FixedSource.hpp"

#include "Particle.hpp"

#include <cstddef>
#include <future>
#include <iostream>
#include <list>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

// FixedSource

//// public

FixedSource::FixedSource(const pugi::xml_node& root)
    : Driver{root}, source{root.child("problemtype").child("fixedsource")} {};

Estimator FixedSource::Solve() {
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
  return estimators;
}

//// private

Estimator FixedSource::StartWorker() {
  Estimator worker_estimator;
  while (const auto range = chunk_giver.Next()) {
    for (auto h = range->first; h < range->second; h++) {
      // a single integer `h` uniquely determines the history of a particle
      // avoid zero seed with h + 1
      auto p{source.Sample(h + 1)};
      // we choose a list because list::splice is constant time
      std::list<Particle> bank(1, p);
      while (!bank.empty()) {
        auto result = bank.back().Transport(world);
        bank.pop_back();
        worker_estimator += result.estimator;
        bank.splice(bank.begin(), result.banked);
      }
    }
  }
  return worker_estimator;
}
