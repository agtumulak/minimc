#include "KEigenvalue.hpp"

#include "Bank.hpp"
#include "Source.hpp"
#include "TransportMethod.hpp"
#include "pugixml.hpp"

#include <future>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

// KEigenvalue

//// public

KEigenvalue::KEigenvalue(const pugi::xml_node& root)
    : Driver{root}, last_inactive{root.child("problemtype")
                                      .child("keigenvalue")
                                      .attribute("inactive")
                                      .as_uint()},
      last_active{
          last_inactive + root.child("problemtype")
                              .child("keigenvalue")
                              .attribute("active")
                              .as_uint()} {
  Source source{
      root.child("problemtype").child("keigenvalue").child("initialsource")};
  for (RNG::result_type s = 1; s <= batchsize; s++) {
    source_bank.push_back(source.Sample(s));
  }
}

EstimatorSet KEigenvalue::Solve() {
  EstimatorSet solver_estimators = init_estimator_set;
  // Perform inactive cycles
  for (Cycle c = 0; c < last_inactive; c++) {
    std::cout << "===== Cycle " << std::to_string(c) << " =====" << std::endl;
    std::cout << "source bank: " << std::to_string(source_bank.size())
              << std::endl;
    std::vector<std::future<std::tuple<Bank, EstimatorSet>>>
        thread_results;
    // start all the threads
    for (size_t i = 0; i < threads; i++) {
      thread_results.push_back(std::async(&KEigenvalue::StartWorker, this));
    }
    // collect results from each thread
    Bank fission_bank;
    for (auto& thread_result : thread_results) {
      auto [worker_bank, worker_estimators] = thread_result.get();
      // TODO: Fix
      // fission_bank += worker_bank;

    }
    // TODO: Fix
    // // swap fission and source bank
    // std::swap(cycle_outcome.banked, source_bank);
  }
}

//// private

std::tuple<Bank, EstimatorSet> KEigenvalue::StartWorker() {
  EstimatorSet worker_estimators = init_estimator_set;
  Bank worker_bank;
  while (auto p = NextParticle()) {
    worker_bank +=
        transport_method->Transport(p.value(), worker_estimators, world);
  }
  return {worker_bank, worker_estimators};
}

std::optional<Particle> KEigenvalue::NextParticle() noexcept {
  std::lock_guard guard(source_bank_mutex);
  if (!source_bank.empty()) {
    auto source_particle = source_bank.back();
    source_bank.pop_back();
    return source_particle;
  }
  return std::nullopt;
}
