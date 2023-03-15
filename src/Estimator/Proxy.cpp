#include "Estimator/Proxy.hpp"

#include "Estimator/Estimator.hpp"
#include "Estimator/Visitor.hpp"
#include "Perturbation/Sensitivity/Proxy/Visitor.hpp"
#include "Perturbation/Sensitivity/Sensitivity.hpp"

#include <type_traits>

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
    // update pending scores
    // https://stackoverflow.com/a/12965621/5101335
    pending_scores[bin] += score;
    // update all sensitivty proxies
    for (auto& sensitivity_proxy : sensitivity_proxies) {
      sensitivity_proxy->Visit(
          *estimator.GetSensitivityProxyVisitor(visitor.particle, bin, score));
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
