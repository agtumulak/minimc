#include "IndirectEffect.hpp"

// IndirectEffect::Visitor

//// public

IndirectEffect::Visitor::~Visitor() noexcept {}

// IndirectEffect

//// public

IndirectEffect::~IndirectEffect() noexcept {}

// TotalCrossSectionPerturbationIndirectEffect

//// public

TotalCrossSectionPerturbationIndirectEffect::
    TotalCrossSectionPerturbationIndirectEffect(
        const TotalCrossSectionPerturbation& perturbation) noexcept
    : perturbation{perturbation} {}

void TotalCrossSectionPerturbationIndirectEffect::Visit(
    const Visitor& visitor) noexcept {
  visitor.Visit(*this);
}

std::unique_ptr<IndirectEffect>
TotalCrossSectionPerturbationIndirectEffect::Clone() const noexcept {
  return std::make_unique<TotalCrossSectionPerturbationIndirectEffect>(
      perturbation);
}
