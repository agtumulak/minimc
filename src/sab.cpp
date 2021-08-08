#include "BasicTypes.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "ThermalScattering.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

/// @brief Samples S(a,b)
/// @details Command line arguments in order
///          (size_t) number of samples
///          (double) incident neutron energy in MeV
///          (double) target temperature in K
///          (string) path to XML file containing only a `tsl` node
///          (string) path to beta boundaries
///          (string) path to alpha boundaries
int main(int argc, char* argv[]) {
  std::cout << "S(a,b) sampler..." << std::endl;
  if (argc != 7) {
    throw std::runtime_error("Expected 6 arguments!");
  }
  const size_t samples{std::stoul(argv[1])};
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
  std::vector<size_t> counts(
      (b_boundaries.size() - 1) * (a_boundaries.size() - 1), 0);

  // samples in one percent
  const size_t one_percent = samples * 0.01;
  for (RNG::result_type seed = 1; seed <= samples; seed++) {
    if (seed % one_percent == 0) {
      std::cout << "\r" << seed / one_percent << "%..." << std::flush;
    }
    Particle p{Point{}, Direction{1, 0, 0}, E, Particle::Type::neutron, seed};
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
    counts[(a_bin - 1) * b_boundaries.size() + (b_bin - 1)]++;
  }
  std::cout << std::endl;
  const std::filesystem::path output_path{"./sab_values.out"};
  std::ofstream outfile{output_path};
  for (auto count : counts) {
    outfile << count << std::endl;
  }
  std::cout << "output written to: " << output_path << std::endl;
}
