#include "Estimator/Estimator.hpp"

#include "Estimator/Visitor.hpp"
#include "Perturbation/Sensitivity/Proxy/Proxy.hpp"
#include "Perturbation/Sensitivity/Proxy/Visitor.hpp"
#include "Perturbation/Sensitivity/Sensitivity.hpp"
#include "World.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <sstream>
#include <string>

using namespace Estimator;

// Interface

//// public

std::unique_ptr<Interface> Interface::Create(
    const pugi::xml_node& estimator_node, const World& world,
    const std::vector<std::unique_ptr<const Perturbation::Interface>>&
        perturbations) {
  const std::string estimator_type = estimator_node.name();
  if (estimator_type == "current") {
    return std::make_unique<Current>(estimator_node, world, perturbations);
  }
  else {
    assert(false); // this should have been caught by the validator
    return {};
  }
}

Interface::~Interface() noexcept {}

void Interface::AddScore(const BinIndex i, const Real s) noexcept {
  sums[i] += s;
  sum_squares[i] += s * s;
}

Real Interface::GetScore(
    const BinIndex i, const Real total_weight) const noexcept {
  return sums.at(i) / total_weight;
}

std::string Interface::to_string(const Real total_weight) const noexcept {
  std::string result;
  // add header
  result += name + "\n" + std::string(name.size(), '=') + "\n\n";
  // add bins
  result += bins->to_string();
  // add mean
  std::stringstream sstream;
  sstream << "mean\n----\n";
  sstream << std::scientific;
  for (const auto& sum : sums) { // wow so clean
    sstream << sum / total_weight << ", ";
  }
  sstream << "\n\n";
  // add standard deviation
  sstream << "std dev\n-------\n";
  for (BinIndex i = 0; i < sums.size(); i++) { // gross
    const auto& sum = sums[i];
    const auto& sum_square = sum_squares[i];
    sstream << std::sqrt(sum_square - sum * sum / total_weight) / total_weight
            << ", ";
  }
  sstream << "\n\n";
  result += sstream.str();
  // add sensitivities
  for (const auto& sensitivity : sensitivities) {
    result += sensitivity->to_string(total_weight);
  }
  return result;
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
  // add perturbation sensitivites
  for (auto sensitivity_it = sensitivities.cbegin(),
            other_it = other.sensitivities.cbegin();
       sensitivity_it != sensitivities.cend(); sensitivity_it++, other_it++) {
    Perturbation::Sensitivity::Interface& sensitivity = **sensitivity_it;
    const Perturbation::Sensitivity::Interface& other = **other_it;
    sensitivity += other;
  }
  return *this;
}

//// protected

Interface::Interface(
    const pugi::xml_node& estimator_node,
    const std::vector<std::unique_ptr<const Perturbation::Interface>>&
        perturbations) noexcept
    : name{estimator_node.attribute("name").as_string()},
      bins{std::make_shared<ParticleBins>(estimator_node.child("bins"))},
      // IIFE
      sensitivities{[this, &estimator_node, &perturbations]() {
        std::vector<std::unique_ptr<Perturbation::Sensitivity::Interface>>
            result;
        for (const auto& sensitivity_node :
             estimator_node.child("perturbations")) {
          result.push_back(Perturbation::Sensitivity::Interface::Create(
              sensitivity_node, perturbations, *this));
        }
        return result;
      }()} {}

Interface::Interface(const Interface& other) noexcept
    : name{other.name}, bins{other.bins}, sensitivities{[&other]() {
        std::vector<std::unique_ptr<Perturbation::Sensitivity::Interface>>
            result;
        for (const auto& sensitivity : other.sensitivities) {
          result.emplace_back(sensitivity->Clone());
        }
        return result;
      }()},
      sums{other.sums}, sum_squares{other.sum_squares} {}

// Current

//// public

Current::Current(
    const pugi::xml_node& current_node, const World& world,
    const std::vector<std::unique_ptr<const Perturbation::Interface>>&
        perturbations)
    : Interface{current_node, perturbations},
      surface{world.FindSurfaceByName(
          current_node.attribute("surface").as_string())} {}

std::unique_ptr<Interface> Current::Clone() const noexcept {
  return std::make_unique<Current>(*this);
}

Score Current::Visit(const Visitor& visitor) const noexcept {
  return visitor.Visit(*this);
}

std::unique_ptr<const Perturbation::Sensitivity::Proxy::Visitor>
Current::GetSensitivityProxyVisitor(
    const Particle& p, const BinIndex i, const Score s,
    const bool add_to_existing) const noexcept {
  class Visitor : public Perturbation::Sensitivity::Proxy::Visitor {
  public:
    Visitor(
        const Particle& p, const BinIndex i, const Score s,
        const bool add_to_existing) noexcept
        : Perturbation::Sensitivity::Proxy::Visitor{
              p, i, s, add_to_existing} {};
    void Visit(Perturbation::Sensitivity::Proxy::TotalCrossSection& proxy)
        const noexcept final {
      // perturbing the total cross section only has an indirect effect
      const auto& indirect_effect = proxy.GetIndirectEfects().front();
      proxy.pending_scores[index] += indirect_effect * score;
    }
  };
  return std::make_unique<Visitor>(p, i, s, add_to_existing);
}
