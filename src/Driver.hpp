#pragma once

#include "BasicTypes.hpp"
#include "World.hpp"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <vector>

namespace Estimator {
class Interface;
class Proxy;
} // namespace Estimator
namespace Perturbation {
namespace IndirectEffect {
class Visitor;
}
class Interface;
} // namespace Perturbation
namespace pugi {
class xml_node;
}
class Particle;
class StreamDelegate;

/// @brief A Driver owns all the data needed to perform radiation transport
class Driver {
public:
  /// @brief Factory method to create new Driver from XML document
  /// @returns A `std::unique_ptr` to the constructed Driver (C++ Core
  ///          Guidelines R.30)
  static std::unique_ptr<Driver> Create(
      const std::filesystem::path& xml_filepath,
      const std::filesystem::path& output_filepath);
  /// @brief Creates objects necessary for a radiation transport
  Driver(
      const pugi::xml_node& root, const std::filesystem::path& output_filepath);
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Driver() noexcept;
  /// @brief Solves the problem, update estimators, and print results
  virtual const std::vector<std::unique_ptr<Estimator::Interface>>&
  Solve() const = 0;

protected:
  /// @brief Kernel of the Monte Carlo transport loop
  /// @details Samples new Particle states until it dies, scoring any
  ///          Estimator::Proxy objects along the way. A vector of
  ///          Estimator::Proxy is passed as an argument since the history may
  ///          not be fully sampled when this function returns.
  void Transport(Particle& p, std::vector<Estimator::Proxy>& estimator_proxies)
      const noexcept;
  /// @brief Global, read-only description of geometric and material properties
  const World world;
  /// @brief Perturbations whose effect on one or more estimators is to be
  ///        estimated
  const std::vector<std::unique_ptr<const Perturbation::Interface>>
      perturbations;
  /// @brief Used to initialize each thread's own estimator as well as
  ///        accumulate results from each thread
  const std::vector<std::unique_ptr<Estimator::Interface>> estimators;
  /// @brief Composition over inheritance delegate for sampling next position
  const std::unique_ptr<const StreamDelegate> stream_delegate;
  /// @brief Path to save results
  const std::filesystem::path output_filepath;
  /// @brief Total histories for fixed-source; cycle weight for k-eigenvalue
  ///        (C++ Core Guidelines C.131)
  const size_t total_weight;
  /// @brief Number of threads dedicated to particle transport
  const size_t threads;
  /// @brief Histories are assigned a seed in [seed, seed + batchsize)
  const RNG::result_type seed;

private:
  // Interface for updating each indirect effect after colliding at the current
  // position.
  std::unique_ptr<const Perturbation::IndirectEffect::Visitor>
  GetCollideWithinCellIndirectEffectVisitor(const Particle& p) const noexcept;
};
