#include "FixedSource.hpp"

#include "Particle.hpp"
#include "TransportMethod.hpp"
#include "pugixml.hpp"

#include <cstddef>
#include <future>
#include <iostream>
#include <list>
#include <memory>
#include <string>
#include <vector>

// FixedSource

//// public

FixedSource::FixedSource(const pugi::xml_node& root)
    : Driver{root}, source{root.child("problemtype").child("fixedsource")} {};

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
  std::cout << '\r' << "Done!          " << std::endl;
  return solver_estimator_set;
}

//// private

EstimatorSet FixedSource::StartWorker() {
  // each thread gets its own EstimatorSet
  EstimatorSet worker_estimator_set = init_estimator_set;
  // keep sampling histories until batchsize histories have been sampled
  while (true) {
    // we _could_ return a newly constructed EstimatorSet from each call to
    // Transport() to keep things more functional, but I suspect it will be
    // super slow, so we use a lightweight EstimatorSet::Proxy
    auto scoring_proxy = worker_estimator_set.GetProxy();
    // atomically update thread-local count of histories started
    auto elapsed = histories_elapsed++;
    if (elapsed >= batchsize) {
      break;
    }
    if (const auto interval_period = batchsize / 10000;
        interval_period != 0 && elapsed % interval_period == 0) {
      std::cout << '\r' << static_cast<double>(elapsed) / batchsize * 100
                << "% complete...";
    }
    std::list<Particle> bank;
    // use the remaining number of histories run for the seed
    bank.emplace_back(source.Sample(seed + elapsed));
    // beginning with the initial Particle, complete the history
    while (!bank.empty()) {
      transport_method->Transport(bank.back(), scoring_proxy, world);
      // new Particle objects are added to front
      bank.back().MoveSecondariesToFrontOf(bank);
      // most recently processed Particle is removed from back
      bank.pop_back();
    }
    // add the scores stored by the proxy to worker_estimator_set
    scoring_proxy.CommitHistory();
  }
  return worker_estimator_set;
}
