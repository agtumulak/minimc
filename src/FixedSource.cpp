#include "FixedSource.hpp"

#include "Bank.hpp"
#include "Particle.hpp"
#include "pugixml.hpp"

#include <future>
#include <iostream>
#include <string>
#include <vector>

// FixedSource

//// public

FixedSource::FixedSource(const pugi::xml_node& root)
    : Driver{root}, source{
                        root.child("problemtype").child("fixedsource"), world,
                        perturbations} {};

EstimatorSet FixedSource::Solve() {
  std::vector<std::future<EstimatorSet>> worker_estimator_sets;
  std::cout << "Spawning " << std::to_string(threads) << " threads working on "
            << std::to_string(batchsize) << " histories..." << std::endl;
  for (size_t i = 0; i < threads; i++) {
    worker_estimator_sets.push_back(
        std::async(&FixedSource::StartWorker, this));
  }
  EstimatorSet solver_estimator_set = init_estimator_set;
  for (auto& worker_estimator_set : worker_estimator_sets) {
    solver_estimator_set += worker_estimator_set.get();
  }
  std::cout << '\r' << "Done!                   " << std::endl;
  return solver_estimator_set;
}

//// private

EstimatorSet FixedSource::StartWorker() {
  // each thread gets its own EstimatorSet
  EstimatorSet worker_estimator_set = init_estimator_set;
  // keep sampling histories until batchsize histories have been sampled
  while (true) {
    // atomically update thread-local count of histories started
    auto elapsed = histories_elapsed++;
    if (elapsed >= batchsize) {
      break;
    }
    if (const auto interval_period = batchsize / 10000;
        interval_period != 0 && elapsed % interval_period == 0) {
      std::cout << '\r' << static_cast<double>(elapsed) / batchsize * 100
                << "% complete...                ";
    }
    // history will be considered fully sampled when bank is empty...
    auto bank = source.Sample(seed + elapsed, seed + elapsed + 1);
    // ...until then we accumulate scores in EstimatorProxy objects
    auto estimator_proxies = worker_estimator_set.CreateProxies();
    // sample the full history
    while (!bank.empty()) {
      // transport the Particle
      Transport(bank.back(), estimator_proxies);
      // new Particle objects are added to front
      bank.back().MoveSecondariesTo(bank);
      // most recently processed Particle is removed from back
      bank.pop_back();
    }
    // history is fully sampled
    for (const auto& estimator_proxy : estimator_proxies) {
      estimator_proxy.CommitHistory();
    }
  }
  return worker_estimator_set;
}
