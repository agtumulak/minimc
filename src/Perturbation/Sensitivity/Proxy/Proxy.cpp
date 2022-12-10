#include "Perturbation/Sensitivity/Proxy/Proxy.hpp"

#include "Particle.hpp"
#include "Perturbation/IndirectEffect/IndirectEffect.hpp"
#include "Perturbation/Sensitivity/Proxy/Visitor.hpp"
#include "Perturbation/Sensitivity/Sensitivity.hpp"

#include <algorithm>
#include <cassert>

using namespace Perturbation::Sensitivity::Proxy;

// Interface

//// public

Interface::Interface(Sensitivity::Interface& sensitivity) noexcept
    : sensitivity{sensitivity} {}

Interface::~Interface() noexcept {}

void Interface::SetIndirectEffects(const Particle& p) noexcept {
  const auto& indirect_effect_it = std::find_if(
      p.indirect_effects.cbegin(), p.indirect_effects.cend(),
      [this](const auto& indirect_effect_ptr) {
        return &indirect_effect_ptr->perturbation == &sensitivity.perturbation;
      });
  if (indirect_effect_it == p.indirect_effects.cend()) {
    assert(false); // particle indirect effects should contain all perturbations
  }
  indirect_effect = *indirect_effect_it;
}

const std::vector<Real>& Interface::GetIndirectEfects() const noexcept {
  return indirect_effect->indirect_effects;
}

// TotalCrossSection

//// public

TotalCrossSection::TotalCrossSection(
    Sensitivity::Interface& sensitivity) noexcept
    : Interface{sensitivity} {}

void TotalCrossSection::Visit(const Visitor& visitor) noexcept {
  visitor.Visit(*this);
}

void TotalCrossSection::CommitHistory() const noexcept {
  for (const auto& [index, score] : pending_scores) {
    sensitivity.AddScore(index, score);
  }
}
