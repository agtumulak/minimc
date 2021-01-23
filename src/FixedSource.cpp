#include "FixedSource.hpp"

#include "History.hpp"

#include <future>
#include <iostream>
#include <numeric>
#include <sstream>

// ChunkGiver

//// public

ChunkGiver::ChunkGiver(RNG::result_type last, size_t chunksize)
    : last{last}, chunksize{chunksize} {}

std::optional<std::pair<RNG::result_type, RNG::result_type>>
ChunkGiver::Next() {
  const std::lock_guard<std::mutex> lock(m);
  next_begin = next_end;
  if (next_begin > last) {
    return std::nullopt;
  }
  next_end = next_begin + chunksize;
  if (next_end > last + 1) {
    return std::make_pair(next_begin, last + 1);
  }
  return std::make_pair(next_begin, next_end);
}

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
    if (!range) {
      break;
    }
    for (auto h = range->first; h < range->second; h++) {
      auto source_particle{Sample(h)};
      worker_estimator += History{source_particle, world}.Transport().estimator;
    }
  }
  return worker_estimator;
}

//// private

Particle FixedSource::Sample(RNG::result_type history) const noexcept {
  RNG rng{history};
  auto p = source.Sample(rng);
  p.seed = rng();
  p.SetCell(world.FindCellContaining(p.GetPosition()));
  return p;
}
