#include "BasicTypes.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "ThermalScattering.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

std::atomic<size_t> histories_elapsed{0};

/// @brief Worker
/// @param sab Thermal scattering data class
/// @param E Incident neutron energy in MeV
/// @param T Target nuclide temperature in K
/// @param max_history History at which workers should stop working
/// @param b_boundaries Bin boundaries for beta
/// @param a_boundaries Bin boundaries for alpha
/// @param progress_interval Number of histories that must pass before progress
///                          bar is updated
std::vector<size_t> Worker(
    const ThermalScattering& sab, const ContinuousEnergy E, const Temperature T,
    const std::vector<double>& b_boundaries, const std::vector<double>& a_boundaries,
    const size_t max_history, const size_t progress_interval) {
  std::vector<size_t> local_counts(
      (b_boundaries.size() - 1) * (a_boundaries.size() - 1), 0);
  while (true) {
    const RNG::result_type history = ++histories_elapsed;
    if (history > max_history) {
      break;
    }
    if (history % progress_interval == 0) {
      std::cout << "\r" << history / progress_interval << "%..." << std::flush;
    }
    Particle p{
        Point{}, Direction{1, 0, 0}, E, Particle::Type::neutron, history};
    const auto beta = sab.SampleBeta(p, E, T);
    const auto alpha = sab.SampleAlpha(p, beta, E, T);
    const size_t b_bin = std::distance(
        b_boundaries.cbegin(),
        std::upper_bound(b_boundaries.cbegin(), b_boundaries.cend(), beta));
    const size_t a_bin = std::distance(
        a_boundaries.cbegin(),
        std::upper_bound(a_boundaries.cbegin(), a_boundaries.cend(), alpha));
    if (b_bin == 0 || a_bin == 0) {
      continue;
    }
    local_counts[(a_bin - 1) * (b_boundaries.size() - 1) + (b_bin - 1)]++;
  }
  return local_counts;
}

/// @brief Samples S(a,b)
/// @details Command line arguments in order
///          (size_t) number of samples
///          (double) incident neutron energy in MeV
///          (double) target temperature in K
///          (string) path to XML file containing only a `tsl` node
///          (string) path to beta boundaries
///          (string) path to alpha boundaries
///          (size_t) number of threads
///          (string) path to output file
int main(int argc, char* argv[]) {
  std::cout << "S(a,b) sampler..." << std::endl;
  if (argc != 9) {
    throw std::runtime_error("Expected 8 arguments!");
  }
  const size_t samples(std::stod(argv[1]));
  const ContinuousEnergy E{std::stod(argv[2])};
  const Temperature T{std::stod(argv[3])};
  // load S(a,b) data
  pugi::xml_document doc;
  doc.load_file(argv[4]);
  const ThermalScattering sab{
      doc.root().select_node("/nuclide/neutron/scatter/tsl").node()};
  // load beta and alpha boundaries
  double value;
  std::vector<double> b_boundaries;
  std::ifstream b_bounds_file{argv[5]};
  while (b_bounds_file >> value) {
    b_boundaries.push_back(value);
  }
  std::vector<double> a_boundaries;
  std::ifstream a_bounds_file{argv[6]};
  while (a_bounds_file >> value) {
    a_boundaries.push_back(value);
  }
  const size_t threads{std::stoul(argv[7])};
  const size_t progress_interval = samples * 0.01;
  const std::filesystem::path output_path{argv[8]};
  // spawn workers
  std::vector<std::future<std::vector<size_t>>> results;
  for (size_t i = 0; i < threads; i++) {
    results.push_back(std::async(
        &Worker, sab, E, T, b_boundaries, a_boundaries, samples,
        progress_interval));
  }
  // add up results
  const auto counts = std::reduce(
      results.begin(), results.end(),
      std::vector<size_t>(
          (b_boundaries.size() - 1) * (a_boundaries.size() - 1), 0),
      [](auto& accumulated, auto& future) {
        std::transform(
            accumulated.cbegin(), accumulated.cend(), future.get().cbegin(),
            accumulated.begin(), std::plus<size_t>{});
        return accumulated;
      });
  std::cout << std::endl;
  std::ofstream outfile{output_path};
  for (auto count : counts) {
    outfile << count << std::endl;
  }
  std::cout << "output written to: " << output_path << std::endl;
}
