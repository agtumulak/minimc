#pragma once

#include "BasicTypes.hpp"
#include "Estimator.hpp"
#include "Perturbation.hpp"
#include "World.hpp"

#include <cstddef>
#include <filesystem>
#include <memory>

namespace pugi {
class xml_node;
}
class TransportMethod;

/// @brief A Driver owns all the data needed to perform radiation transport
class Driver {
public:
  /// @brief Factory method to create new Driver from XML document
  /// @returns A `std::unique_ptr` to the constructed Driver (C++ Core
  ///          Guidelines R.30)
  static std::unique_ptr<Driver>
  Create(const std::filesystem::path& xml_filepath);
  /// @brief Creates objects necessary for a radiation transport
  Driver(const pugi::xml_node& root);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Driver() noexcept;
  /// @brief Solves the problem
  virtual EstimatorSet Solve() = 0;

protected:
  /// @brief Global, read-only description of geometric and material properties
  const World world;
  /// @brief All perturbations whose effect on one or more Estimator objects is
  ///        being estimated
  const PerturbationSet perturbations;
  /// @brief Total histories for fixed-source; cycle weight for k-eigenvalue
  const RNG::result_type batchsize;
  /// @brief Initial value of an EstimatorSet used to initialize workers
  const EstimatorSet init_estimator_set;
  /// @brief Number of threads dedicated to particle transport
  const size_t threads;
  /// @brief Histories are assigned a seed in [seed, seed + batchsize)
  const RNG::result_type seed;
};
