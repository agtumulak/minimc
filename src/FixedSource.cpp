#include "FixedSource.hpp"

#include "Particle.hpp"

#include <cstddef>
#include <future>
#include <iostream>
#include <list>
#include <numeric>
#include <string>
#include <vector>

// FixedSource

//// public

FixedSource::FixedSource(const pugi::xml_node& root)
    : Driver{root}, source{root.child("problemtype").child("fixedsource")} {};

Estimator FixedSource::Solve() {
  std::vector<std::future<Estimator>> results;
  std::cout << "Spawning " << std::to_string(threads) << " threads working on "
            << std::to_string(batchsize) << " histories..." << std::endl;
  for (size_t i = 0; i < threads; i++) {
    results.push_back(std::async(&FixedSource::StartWorker, this));
  }
  estimators = std::reduce(
      results.begin(), results.end(), Estimator{},
      [](auto& accumulated, auto& future) {
        return accumulated += future.get();
      });
  estimators.Normalize(batchsize);
  return estimators;
}

//// private

Estimator FixedSource::StartWorker() {
  Estimator worker_estimator;
  while (true) {
    // atomically update thread-local count of histories started
    auto elapsed = histories_elapsed++;
    if (elapsed >= batchsize) {
      break;
    }
    // Use the remaining number of histories remaining as the seed. This is
    // guaranteed to be greater than one.
    auto p{source.Sample(batchsize - seed)};
    // we choose a list because list::splice is constant time
    std::list<Particle> bank(1, p);
    while (!bank.empty()) {
      auto result = bank.back().Transport(world);
      bank.pop_back();
      worker_estimator += result.estimator;
      bank.splice(bank.begin(), result.banked);
    }
  }
  return worker_estimator;
}
