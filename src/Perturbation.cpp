#include "Perturbation.hpp"

#include <stdexcept>
#include <string>

Perturbation ToPerturbation(const std::string& name) {
  if (name == "total") {
    return Perturbation::total;
  }
  else if (name == "capture") {
    return Perturbation::capture;
  }
  else if (name == "scatter") {
    return Perturbation::scatter;
  }
  else {
    throw std::runtime_error("Unrecognized perturbation name: " + name);
  };
}
