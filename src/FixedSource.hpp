#pragma once

#include "BasicTypes.hpp"
#include "Driver.hpp"
#include "Particle.hpp"
#include "pugixml.hpp"

#include <mutex>
#include <optional>
#include <utility>

/// @brief Splits histories into fixed-size chunks.
/// @details The last chunk may have fewer values than `chunksize`. Each call
///          to Next() will return a distinct range in a thread-safe manner.
class ChunkGiver {
public:
  /// @brief Constructs a ChunkGiver with a last history and fixed chunk size
  ChunkGiver(RNG::result_type last, size_t chunksize);
  /// @brief Disallow copy/move constructor:
  ///        https://www.stroustrup.com/C++11FAQ.html#default2
  ChunkGiver(const ChunkGiver&) = delete;
  /// @brief Disallow copy/move assignment operator:
  ///        https://www.stroustrup.com/C++11FAQ.html#default2
  ChunkGiver& operator=(const ChunkGiver&) = delete;
  /// @brief Returns the next chunk. Returns std::nullopt if there are zero
  ///        chunks left to give.
  /// @return A pair of `[start, end)` history numbers.
  std::optional<std::pair<RNG::result_type, RNG::result_type>> Next();

private:
  std::mutex m;
  const RNG::result_type last;
  const size_t chunksize;
  RNG::result_type next_begin{1};
  RNG::result_type next_end{1};
};

/// @brief Creates and executes a fixed-source calculation
class FixedSource : public Driver {
public:
  /// @brief Creates objects necessary for a fixed source calculation
  /// @param root Root node of existing XML document
  FixedSource(const pugi::xml_node& root);
  /// @brief Spawn workers to work on chunks of history
  void Solve() override;
  /// @brief Function exiected by a worker on a single thread
  Estimator StartWorker();

private:
  // Helper function to create a set of Particle::Type which will spawn
  // TODO: Support more than one Particle::Type
  static Particle::Type CreateParticleType(const pugi::xml_node& root) noexcept;
  // Helper function to set default initial energy of spawned Particle
  static Energy CreateDefaultEnergy(const pugi::xml_node& root) noexcept;
  // In a fixed-source calculation a single integer (`history`) uniquely
  // specifies the history of a particle.
  Particle Sample(RNG::result_type history) const noexcept;

  ChunkGiver chunk_giver{batchsize, chunksize};
  const Particle::Type particle_type;
  const Energy initial_energy;
};
