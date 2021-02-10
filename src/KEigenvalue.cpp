#include "KEigenvalue.hpp"

#include "Constants.hpp"
#include "Parallel.hpp"

#include <algorithm>
#include <cstdint>
#include <future>
#include <iterator>
#include <map>
#include <numeric>
#include <random>

// KEigenvalue

//// public

KEigenvalue::KEigenvalue(const pugi::xml_node& root)
    : Driver{root}, last_inactive{root.child("problemtype")
                                      .child("keigenvalue")
                                      .attribute("inactive")
                                      .as_uint()},
      last_active{root.child("problemtype")
                      .child("keigenvalue")
                      .attribute("active")
                      .as_uint()} {
  RNG rng;
  Source source{
      root.child("problemtype").child("keigenvalue").child("initialsource")};
  for (RNG::result_type s = 1; s <= batchsize; s++) {
    rng.seed(s * constants::seed_stride);
    source_bank.push_back(source.Sample(rng));
  }
  k = 1.0;
}

void KEigenvalue::Solve() {
  RNG rng{};
  // Perform inactive cycles
  for (Cycle c = 0; c < last_inactive; c++) {
    Estimator cycle_estimator;
    std::cout << "===== Cycle " << std::to_string(c) << " =====" << std::endl;
    std::cout << "source bank: " << std::to_string(source_bank.size())
              << std::endl;
    std::vector<std::future<std::map<size_t, TransportOutcome>>> thread_results;
    chunk_giver.Reset(source_bank.size());
    for (size_t i = 0; i < threads; i++) {
      thread_results.push_back(std::async(&KEigenvalue::StartWorker, this));
    }
    std::map<size_t, TransportOutcome> sorted_chunks;
    for (auto& thread_result : thread_results) {
      sorted_chunks.merge(thread_result.get());
    }
    std::vector<Particle> sorted_fission_bank;
    for (const auto& chunk_result : sorted_chunks) {
      const auto& [index, outcome] = chunk_result;
      std::move(
          outcome.banked.begin(), outcome.banked.end(),
          std::back_insert_iterator(sorted_fission_bank));
      cycle_estimator += outcome.estimator;
    }
    auto fission_weight{cycle_estimator.at(Estimator::Event::implicit_fission)};
    std::cout << "fission weight: " << fission_weight << std::endl;
    k = k * fission_weight / batchsize;
    std::cout << "k: " << std::to_string(k) << std::endl;
    size_t new_bank_size(batchsize / k + std::uniform_real_distribution{}(rng));
    source_bank.clear();
    std::uniform_int_distribution<size_t> sample_index{
        0, sorted_fission_bank.size() - 1};
    for (size_t i = 0; i < new_bank_size; i++) {
      source_bank.push_back(sorted_fission_bank.at(sample_index(rng)));
    }
  }
}

std::map<size_t, TransportOutcome> KEigenvalue::StartWorker() {
  std::map<size_t, TransportOutcome> chunk_outputs;
  TransportOutcome chunk_output;
  while (const auto range = chunk_giver.Next()) {
    Bank source;
    std::move(
        std::next(source_bank.begin(), range->first),
        std::next(source_bank.begin(), range->second),
        std::back_insert_iterator(source));
    chunk_outputs[range->first] = TransportAndBank(source, world);
  }
  return chunk_outputs;
}
