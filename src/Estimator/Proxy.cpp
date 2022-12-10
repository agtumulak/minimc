#include "Estimator/Proxy.hpp"

#include "Estimator/Visitor.hpp"
#include "Perturbation/Sensitivity/Proxy/Visitor.hpp"
#include "Perturbation/Sensitivity/Sensitivity.hpp"

using namespace Estimator;

// Proxy

//// public

Proxy::Proxy(Interface& estimator) noexcept
    : estimator{estimator},
      // IIFE
      sensitivity_proxies{[&estimator]() noexcept {
        std::vector<
            std::unique_ptr<Perturbation::Sensitivity::Proxy::Interface>>
            result;
        for (const auto& sensitivity : estimator.sensitivities) {
          result.emplace_back(sensitivity->CreateProxy());
        }
        return result;
      }()} {}

void Proxy::SetIndirectEffects(const Particle& p) noexcept {
  for (auto& sensitivity_proxy : sensitivity_proxies) {
    sensitivity_proxy->SetIndirectEffects(p);
  }
}

void Proxy::Visit(const Visitor& visitor) noexcept {
  const auto score = estimator.Visit(visitor);
  if (score != 0) {
    const auto bin = estimator.bins->GetIndex(visitor.particle);
    // check if the bin has already been scored to
    const auto it = pending_scores.find(bin);
    const auto add_to_existing = it != pending_scores.cend();
    // update pending scores
    if (add_to_existing) {
      it->second += score;
    }
    else {
      pending_scores.insert(std::make_pair(bin, score));
    }
    // update all sensitivty proxies
    for (auto& sensitivity_proxy : sensitivity_proxies) {
      sensitivity_proxy->Visit(*estimator.GetSensitivityProxyVisitor(
          visitor.particle, bin, score, add_to_existing));
    }
  }
}

void Proxy::CommitHistory() const noexcept {
  for (const auto& [index, score] : pending_scores) {
    estimator.AddScore(index, score);
  }
  // commit child ScorableProxy objects for Sensitivity objects
  for (auto& sensitivity_proxy : sensitivity_proxies) {
    sensitivity_proxy->CommitHistory();
  }
}
