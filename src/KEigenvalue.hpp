#pragma once

#include "BasicTypes.hpp"
#include "Driver.hpp"
#include "Parallel.hpp"
#include "Particle.hpp"
#include "pugixml.hpp"

#include <cstddef>
#include <list>
#include <map>

/// @brief Creates and executes a k-eigenvalue calculation
class KEigenvalue : public Driver {
public:
  /// @brief Creates objects necessary for a k-eigenvalue calculation
  /// @param root Root node of existing XML document
  KEigenvalue(const pugi::xml_node& root);
  /// @brief Solve sequential fixed-source calculations using a fission bank
  ///        passed between cycles
  void Solve() override;
  /// @brief Function executed by a worker on a single thread
  /// @details Chunks are kept separate as they will have to be sorted in order
  ///          to keep results deterministic
  std::map<size_t, Particle::TransportOutcome> StartWorker();

private:
  using Cycle = size_t;
  std::list<Particle> source_bank{};
  Real k;
  ChunkGiver chunk_giver{batchsize, chunksize};
  const Cycle last_inactive;
  const Cycle last_active;
};
