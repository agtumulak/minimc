#include "FixedSource.hpp"

#include "Bank.hpp"
#include "Estimator/Proxy.hpp"
#include "Particle.hpp"
#include "pugixml.hpp"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <string>
#include <vector>

// FixedSource

//// public

FixedSource::FixedSource(
    const pugi::xml_node& root, const std::filesystem::path& output_filepath)
    : Driver{root, output_filepath},
      source{
          root.child("problemtype").child("fixedsource"), world,
          perturbations} {};

const std::vector<std::unique_ptr<Estimator::Interface>>&
FixedSource::Solve() const {
  // outer vector corresponds to thread, inner vector corresponds to estimator
  std::vector<std::future<std::vector<std::unique_ptr<Estimator::Interface>>>>
      threads_future_estimators;
  std::cout << "Spawning " << std::to_string(threads) << " threads working on "
            << std::to_string(total_weight) << " histories..." << std::endl;
  // Number of histories completed or initiated by all threads. May exceed
  // batchsize since each thread will call it once before ending.
  std::atomic<uint_fast64_t> histories_elapsed{0};
  for (size_t i = 0; i < threads; i++) {
    // https://stackoverflow.com/q/18359864/5101335
    threads_future_estimators.push_back(std::async(
        &FixedSource::StartWorker, this, std::ref(histories_elapsed)));
  }
  // initialize estimators used to accumulate results from each thread
  for (auto& thread_future_estimators : threads_future_estimators) {
    const auto thread_estimators = thread_future_estimators.get();
    // https://stackoverflow.com/a/57048332/5101335
    for (auto main_estimator_it = estimators.cbegin(),
              thread_estimator_it = thread_estimators.cbegin();
         main_estimator_it != estimators.cend();
         main_estimator_it++, thread_estimator_it++) {
      Estimator::Interface& main_estimator = **main_estimator_it;
      const Estimator::Interface& thread_estimator = **thread_estimator_it;
      main_estimator += thread_estimator;
    }
  }
  std::cout << '\r' << "Done!                   " << std::endl;
  // save estimator results to file
  if (!output_filepath.empty()) {
    std::ofstream output_file{output_filepath};
    output_file << total_weight << std::endl << std::endl;
    for (const auto& estimator : estimators) {
      output_file << estimator->to_string(total_weight);
    }
    std::cout << "Output written to "
              << std::filesystem::absolute(output_filepath) << std::endl;
  }
  return estimators;
}

//// private

std::vector<std::unique_ptr<Estimator::Interface>>
FixedSource::StartWorker(std::atomic<uint_fast64_t>& histories_elapsed) const {
  // each thread gets its own copy of Estimator::Interface objects
  std::vector<std::unique_ptr<Estimator::Interface>> worker_estimators;
  for (const auto& estimator : estimators) {
    worker_estimators.emplace_back(estimator->Clone());
  }
  // keep sampling histories until all histories have been sampled
  while (true) {
    // atomically update thread-local count of histories started
    auto elapsed = histories_elapsed++;
    if (elapsed >= total_weight) {
      break;
    }
    if (const auto interval_period = total_weight / 10000;
        interval_period != 0 && elapsed % interval_period == 0) {
      std::cout << '\r' << static_cast<double>(elapsed) / total_weight * 100
                << "% complete...                ";
    }
    // history will be considered fully sampled when bank is empty...
    auto bank = source.Sample(seed + elapsed, seed + elapsed + 1);
    // ...until then we accumulate scores in Estimator::Proxy objects
    std::vector<Estimator::Proxy> estimator_proxies;
    for (const auto& worker_estimator : worker_estimators) {
      estimator_proxies.emplace_back(*worker_estimator);
    }
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
  return worker_estimators;
}
