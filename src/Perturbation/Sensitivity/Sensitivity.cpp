#include "Perturbation/Sensitivity/Sensitivity.hpp"

#include "Estimator/Estimator.hpp"
#include "Perturbation/Sensitivity/Proxy/Proxy.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace Perturbation::Sensitivity;

// Sensitivity

//// public

std::unique_ptr<Interface> Interface::Create(
    const pugi::xml_node& perturbations_sensitivity_node,
    const std::vector<std::unique_ptr<const Perturbation::Interface>>&
        perturbations,
    const Estimator::Interface& estimator) {
  // find matching Perturbation
  const std::string& name =
      perturbations_sensitivity_node.attribute("to").as_string();
  const auto perturbation_it = std::find_if(
      perturbations.cbegin(), perturbations.cend(),
      [&name](const auto& perturbation_ptr) {
        return perturbation_ptr->name == name;
      });
  if (perturbation_it == perturbations.cend()) {
    throw std::runtime_error(
        "Perturbation \"" + name + "\" not found. Must be one of: [" +
        std::accumulate(
            perturbations.cbegin(), perturbations.cend(), std::string{},
            [](const auto& accumulated, const auto& perturbation_ptr) noexcept {
              return accumulated + "\"" + perturbation_ptr->name + "\", ";
            }) +
        "]");
  }
  return (*perturbation_it)->CreateSensitivity(estimator);
}

Interface::Interface(
    const Estimator::Interface& estimator,
    const Perturbation::Interface& perturbation) noexcept
    : name{estimator.name + "::" + perturbation.name},
      perturbation{perturbation}, bins{*estimator.bins} {}

Interface::~Interface() noexcept {}

void Interface::AddScore(const size_t index, const Score s) noexcept {
  sums[index] += s;
  sum_squares[index] += s * s;
}

Score Interface::GetScore(
    const size_t i, const Real total_weight) const noexcept {
  return sums.at(i) / total_weight;
}

Interface& Interface::operator+=(const Interface& other) noexcept {
  // add scores
  std::transform(
      sums.cbegin(), sums.cend(), other.sums.cbegin(), sums.begin(),
      std::plus<Real>());
  // add sum of squares
  std::transform(
      sum_squares.cbegin(), sum_squares.cend(), other.sum_squares.cbegin(),
      sum_squares.begin(), std::plus<Real>());
  return *this;
}

// TotalCrossSection

//// public

TotalCrossSection::TotalCrossSection(
    const Estimator::Interface& estimator,
    const Perturbation::Interface& perturbation) noexcept
    : Interface{estimator, perturbation} {}

std::unique_ptr<Interface> TotalCrossSection::Clone() const noexcept {
  return std::make_unique<TotalCrossSection>(*this);
}

std::unique_ptr<Proxy::Interface> TotalCrossSection::CreateProxy() noexcept {
  return std::make_unique<Perturbation::Sensitivity::Proxy::TotalCrossSection>(
      *this);
}

std::string
TotalCrossSection::to_string(const Real total_weight) const noexcept {
  std::string result;
  // add header
  result += name + "\n" + std::string(name.size(), '=') + "\n\n";
  // add mean
  std::stringstream sstream;
  sstream << "mean\n----\n";
  sstream << std::scientific;
  for (const auto& sum : sums) {
    sstream << sum / total_weight << ", ";
  }
  sstream << "\n\n";
  // add standard deviation
  sstream << "std dev\n-------\n";
  for (BinIndex i = 0; i < sums.size(); i++) {
    const auto& sum = sums[i];
    const auto& sum_square = sum_squares[i];
    sstream << std::sqrt(sum_square - sum * sum / total_weight) / total_weight
            << ", ";
  }
  sstream << "\n\n";
  result += sstream.str();
  return result;
}

// TNSL

//// public

TNSL::TNSL(
    const Estimator::Interface& estimator,
    const Perturbation::Interface& perturbation) noexcept
    : Interface{estimator, perturbation} {}

std::unique_ptr<Interface> TNSL::Clone() const noexcept {
  return std::make_unique<TNSL>(*this);
}

std::unique_ptr<Proxy::Interface> TNSL::CreateProxy() noexcept {
  return std::make_unique<Perturbation::Sensitivity::Proxy::TNSL>(*this);
}

std::string TNSL::to_string(const Real total_weight) const noexcept {
  std::string result;
  // add header
  result += name + "\n" + std::string(name.size(), '=') + "\n\n";
  // add mean
  std::stringstream sstream;
  sstream << "mean\n----\n";
  sstream << std::scientific;
  for (const auto& sum : sums) {
    sstream << sum / total_weight << ", ";
  }
  sstream << "\n\n";
  // add standard deviation
  sstream << "std dev\n-------\n";
  for (BinIndex i = 0; i < sums.size(); i++) {
    const auto& sum = sums[i];
    const auto& sum_square = sum_squares[i];
    sstream << std::sqrt(sum_square - sum * sum / total_weight) / total_weight
            << ", ";
  }
  sstream << "\n\n";
  result += sstream.str();
  return result;
}
