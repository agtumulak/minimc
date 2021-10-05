#include "Bins.hpp"

#include <cassert>
#include <stdexcept>
#include <string>

// Bins

std::ostream& operator<<(std::ostream& os, const Bins& b) noexcept {
  b.Print(os);
  return os;
}

//// public

std::unique_ptr<const Bins>
Bins::Create(const pugi::xml_node& bins_node) noexcept {
  const std::string bins_type = bins_node.name();
  if (bins_type.empty()) {
    return std::make_unique<const NoBins>();
  }
  else if (bins_type == "linspace") {
    return std::make_unique<const LinspaceBins>(bins_node);
  }
  else {
    assert(false); // this should have been caught by the validator
  }
}

Bins::~Bins() noexcept {}

// NoBins

//// public

size_t NoBins::size() const noexcept { return 1; }

size_t NoBins::GetIndex(const Real&) const noexcept { return 0; }

//// private

void NoBins::Print(std::ostream& os) const noexcept {
  os << "none";
}

// LinspaceBins

//// public

LinspaceBins::LinspaceBins(const pugi::xml_node& linspace_node) noexcept
    : n_bins{linspace_node.attribute("bins").as_ullong()},
      lower_bound{linspace_node.attribute("min").as_double()},
      bin_width{
          (linspace_node.attribute("max").as_double() - lower_bound) /
          (n_bins - 1)} {}

size_t LinspaceBins::size() const noexcept { return n_bins; }

size_t LinspaceBins::GetIndex(const Real& v) const {
  const auto bin_index = (v - lower_bound) / bin_width;
  if (bin_index >= n_bins - 1) {
    throw std::out_of_range(
        "Value (" + std::to_string(v) +
        ") is greater than or equal to largest bin boundary {" +
        std::to_string(lower_bound + (n_bins - 1) * bin_width) + ")");
  }
  else if (bin_index < 0) {
    return 0;
  }
  else {
    return size_t(bin_index) + 1;
  }
}

//// private

void LinspaceBins::Print(std::ostream& os) const noexcept {
  for (size_t i = 0; i < n_bins; i++) {
    os << std::to_string(lower_bound + i * bin_width) << " ";
  }
}
