#include "Perturbation/IndirectEffect/IndirectEffect.hpp"

#include "Perturbation/IndirectEffect/Visitor.hpp"
#include "Perturbation/Perturbation.hpp"

using namespace Perturbation::IndirectEffect;

// Interface

//// public

Interface::Interface(const Perturbation::Interface& perturbation) noexcept
    : perturbation{perturbation} {}

Interface::~Interface() noexcept {}

// TotalCrossSection

//// public

TotalCrossSection::TotalCrossSection(
    const Perturbation::TotalCrossSection& perturbation) noexcept
    : Interface{perturbation}, nuclide{perturbation.nuclide} {}

void TotalCrossSection::Visit(const Visitor& visitor) noexcept {
  visitor.Visit(*this);
}

std::unique_ptr<Interface> TotalCrossSection::Clone() const noexcept {
  return std::make_unique<TotalCrossSection>(*this);
}

// TNSL

//// public

TNSL::TNSL(const Perturbation::TNSL& perturbation) noexcept
    : Interface{perturbation}, perturbation{perturbation} {}

void TNSL::Visit(const Visitor& visitor) noexcept { visitor.Visit(*this); }

std::unique_ptr<Interface> TNSL::Clone() const noexcept {
  return std::make_unique<TNSL>(*this);
}
