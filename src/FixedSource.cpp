#include "FixedSource.hpp"

#include "Nuclide.hpp"
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
  return solver_estimator_set.Normalize(batchsize);
}

//// private

EstimatorSet FixedSource::StartWorker() {
  // we _could_ return a newly constructed EstimatorSet from each call to
  // Transport() to keep it purely functional, but I suspect it will be super
  // slow, so we construct it once at the beginning of a thread and pass it
  // as a reference parameter to Transport()
  EstimatorSet worker_estimator_set = init_estimator_set;
  while (true) {
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
    // Model a random walk as sampling a tree
    // Start with a vector of vector of Histories. The 
    std::vector<History> history_bank;
    // ...and initialize it by sampling from the fixed source, using a seed
    // that should be unique to this history
    history_bank.emplace_back(world, source, seed + elapsed);
    do {
      auto history = history_bank.back();
      history_bank.pop_back();
      // complete the history
      transport_method->Transport(history, world);
      // score any Estimators
      worker_estimator_set.Score(history);
      // Push secondaries onto the bank
      CreateSecondaries(history_bank, history);
    } while (!history_bank.empty());
  }
  return worker_estimator_set;
}

void FixedSource::CreateSecondaries(
    std::vector<History> bank, const History& h) {
  for (const auto& state: h){
    if (state.event == State::Event::fission) {
      state.nuclide->
    }
  }
}
