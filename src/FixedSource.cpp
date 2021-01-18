#include "FixedSource.hpp"

#include "Constants.hpp"
#include "XMLDocument.hpp"

#include <future>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <vector>

// FixedSource

//// public

FixedSource::FixedSource(
    Estimator& estimators, const World& world, const Source& source,
    History histories, size_t threads, size_t chunksize)
    : global{estimators}, world{world}, source{source}, histories{histories},
      threads{threads}, chunksize{chunksize} {}

void FixedSource::PoolSolve() {
  std::cout << "Spawning " << std::to_string(threads) << " threads working on "
            << std::to_string(histories) << " histories split into "
            << std::to_string(
                   histories / chunksize + (histories % chunksize != 0))
            << " chunks... " << std::endl;
  std::vector<std::future<Estimator>> results;
  for (size_t i = 1; i <= threads; i++) {
    results.push_back(std::async(&FixedSource::StartWorker, this));
  }
  const auto result = std::reduce(
      results.begin(), results.end(), Estimator{},
      [](auto& accumulated, auto& future) {
        return accumulated += future.get();
      });
  global += result;
  std::cout << global;
}

//// private

Estimator FixedSource::StartWorker() {
  Estimator worker_estimator;
  while (const auto range = history_chunks.Next()) {
    if (!range) {
      break;
    }
    for (History h = range->first; h < range->second; h++) {
      std::minstd_rand rng{h};
      auto p{source.Sample(rng, world)};
      while (p.IsAlive()) {
        Real collision_distance, surfacecross_distance;
        std::shared_ptr<const CSGSurface> nearest_surface;
        collision_distance = p.GetCell().SampleCollisionDistance(rng, p);
        std::tie(nearest_surface, surfacecross_distance) =
            p.GetCell().NearestSurface(p.GetPosition(), p.GetDirection());
        if (collision_distance < surfacecross_distance) {
          p.Stream(collision_distance);
          const auto& nuclide{p.GetCell().material->SampleNuclide(rng, p)};
          switch (nuclide.SampleReaction(rng, p)) {
          case NuclearData::Reaction::capture:
            worker_estimator.at(Estimator::Event::capture)++;
            p.Kill();
            break;
          case NuclearData::Reaction::scatter:
            worker_estimator.at(Estimator::Event::scatter)++;
            nuclide.Scatter(rng, p);
            break;
          }
          worker_estimator.at(Estimator::Event::collision)++;
        }
        else {
          p.Stream(surfacecross_distance + constants::nudge);
          worker_estimator.at(Estimator::Event::surface_crossing)++;
          p.SetCell(world.FindCellContaining(p.GetPosition()));
          if (!p.GetCell().material) {
            p.Kill();
          }
        }
      }
    }
  }
  return worker_estimator;
}

// FixedSourceStandalone

//// public

FixedSourceStandalone::FixedSourceStandalone(const pugi::xml_node& root)
    : world{root}, source{root},
      histories(
          std::stoi(root.child("general").child("histories").child_value())),
      threads{
          std::stoul(root.child("general").child("histories").child_value())},
      chunksize{
          std::stoul(root.child("general").child("chunksize").child_value())} {}

void FixedSourceStandalone::Solve() {
  FixedSource{global, world, source, histories, threads, chunksize}.PoolSolve();
}
