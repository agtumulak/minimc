#include "FixedSource.hpp"

#include "Constants.hpp"

#include <algorithm>
#include <cstddef>
#include <future>
#include <string>
#include <iostream>
#include <iterator>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>

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
    for (auto h = range->first; h < range->second; h++) {
      std::vector<Particle> bank(1, Sample(h));
      while (!bank.empty()) {
        auto result = bank.back().Transport(world);
        bank.pop_back();
        worker_estimator += result.estimator;
        std::move(
            result.banked.begin(), result.banked.end(),
            std::back_insert_iterator(bank));
      }
    }
  }
  return worker_estimator;
}

//// private

Particle FixedSource::Sample(RNG::result_type history) const noexcept {
  // avoid zero seed with +1
  RNG rng{(history + 1) * constants::seed_stride};
  return source.Sample(rng);
}
