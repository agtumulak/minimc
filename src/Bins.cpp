#include "Bins.hpp"

#include "Particle.hpp"
#include "pugixml.hpp"

#include <cassert>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>

// Bins

//// public

std::unique_ptr<const Bins> Bins::Create(const pugi::xml_node& bins_node) {
  // construct concrete type
  const std::string bins_type = bins_node.name();
  if (bins_type.empty()) {
    return std::make_unique<const NoBins>();
  }
  // check that min < max
  if (!(bins_node.attribute("min").as_double() <
        bins_node.attribute("max").as_double())) {
    throw std::runtime_error(
        bins_node.path() + ": condition not satisfied: min < max");
  }
  if (bins_type == "linspace") {
    return std::make_unique<const LinspaceBins>(bins_node);
  }
  else if (bins_type == "logspace"){
    return std::make_unique<const LogspaceBins>(bins_node);
  }
  else {
    assert(false); // this should have been caught by the validator
    return {};
  }
}

Bins::~Bins() noexcept {}

// NoBins

//// public

size_t NoBins::size() const noexcept { return 1; }

size_t NoBins::GetIndex(const Real&) const noexcept { return 0; }

std::string NoBins::to_string() const noexcept { return "none"; }

// LinspaceBins

//// public

LinspaceBins::LinspaceBins(const pugi::xml_node& linspace_node) noexcept
    : n_bins{linspace_node.attribute("bins").as_ullong() + 2},
      lower_bound{linspace_node.attribute("min").as_double()},
      upper_bound{linspace_node.attribute("max").as_double()},
      bin_width{(upper_bound - lower_bound) / (n_bins - 2)} {}

size_t LinspaceBins::size() const noexcept { return n_bins; }

size_t LinspaceBins::GetIndex(const Real& v) const noexcept {
  if (v < lower_bound){
    return 0;
  }
  else if (v >= upper_bound){
    return n_bins - 1;
  }
  else {
    return (v - lower_bound) / bin_width + 1;
  }
}

std::string LinspaceBins::to_string() const noexcept {
  std::stringstream sstream;
  for (size_t i = 0; i < n_bins - 2; i++) {
    sstream << std::scientific << lower_bound + i * bin_width << ", ";
  }
  return sstream.str();
}

// LogspaceBins

//// public

LogspaceBins::LogspaceBins(const pugi::xml_node& logspace_node) noexcept
    : n_bins{logspace_node.attribute("bins").as_ullong() + 2},
      base{logspace_node.attribute("base").as_double(10)},
      log_lower_bound{logspace_node.attribute("min").as_double()},
      log_upper_bound{logspace_node.attribute("max").as_double()},
      log_bin_width{(log_upper_bound - log_lower_bound) / (n_bins - 2)} {}

size_t LogspaceBins::size() const noexcept { return n_bins; }

size_t LogspaceBins::GetIndex(const Real& v) const noexcept {
  const auto log_v = std::log(v) / std::log(base);
  if (log_v < log_lower_bound) {
    return 0;
  }
  else if (log_v >= log_upper_bound) {
    return n_bins - 1;
  }
  else {
    return (log_v - log_lower_bound) / log_bin_width + 1;
  }
}

std::string LogspaceBins::to_string() const noexcept {
  std::stringstream sstream;
  for (size_t i = 0; i < n_bins - 2; i++) {
    sstream << std::scientific
            << std::pow(base, log_lower_bound + i * log_bin_width) << ", ";
  }
  return sstream.str();
}

// ParticleBins

//// public

ParticleBins::ParticleBins(const pugi::xml_node& bins_node) noexcept
    : // IIFE
      direction{[&bins_node]() noexcept {
        const auto& cosine_node = bins_node.child("cosine");
        return cosine_node
                   ? std::optional<Direction>(
                         std::in_place,
                         cosine_node.attribute("u").as_double(),
                         cosine_node.attribute("v").as_double(),
                         cosine_node.attribute("w").as_double())
                   : std::nullopt;
      }()},
      cosine{Bins::Create(bins_node.child("cosine").first_child())},
      energy{Bins::Create(bins_node.child("energy").first_child())},
      strides{ComputeStrides(*cosine, *energy)} {}

size_t ParticleBins::size() const noexcept {
  return cosine->size() * strides.front();
}

size_t ParticleBins::GetIndex(const Particle& p) const noexcept {
  // get cosine bin
  const auto c_i =
      direction ? cosine->GetIndex(direction.value().Dot(p.GetDirection())) : 0;
  // get energy bin
  const auto e_i = energy->GetIndex(std::visit(VisitEnergy(), p.GetEnergy()));
  // get flattened index
  return strides[0] * c_i + e_i;
}

std::string ParticleBins::to_string() const noexcept {
  std::string result;
  result += "cosine\n------\n" + cosine->to_string() + "\n\n";
  result += "energy\n------\n" + energy->to_string() + "\n\n";
  return result;
}
