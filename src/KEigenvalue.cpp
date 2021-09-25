#include "KEigenvalue.hpp"

#include "Source.hpp"

#include <future>
#include <iostream>
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

Estimator KEigenvalue::Solve() {
  Estimator estimator;
  // Perform inactive cycles
  for (Cycle c = 0; c < last_inactive; c++) {
    std::cout << "===== Cycle " << std::to_string(c) << " =====" << std::endl;
    std::cout << "source bank: " << std::to_string(source_bank.size())
              << std::endl;
    std::vector<std::future<TransportMethod::Outcome>> thread_results;
    // start all the threads
    for (size_t i = 0; i < threads; i++) {
      thread_results.push_back(std::async(&KEigenvalue::StartWorker, this));
    }
    // collect results from each thread
    TransportMethod::Outcome cycle_outcome;
    for (auto& thread_result : thread_results) {
      cycle_outcome += thread_result.get();
    }
    // swap fission and source bank
    std::swap(cycle_outcome.banked, source_bank);
  }
  return estimator;
}

//// private

TransportMethod::Outcome KEigenvalue::StartWorker() {
  TransportMethod::Outcome worker_outcome;
  while (auto p = NextParticle()) {
    worker_outcome += transport_method->Transport(p.value(), world);
  }
  return worker_outcome;
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
