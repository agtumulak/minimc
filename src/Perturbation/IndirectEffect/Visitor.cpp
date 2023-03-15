#include "Perturbation/IndirectEffect/Visitor.hpp"

using namespace Perturbation::IndirectEffect;

// Visitor

//// public

Visitor::~Visitor() noexcept {}

void Visitor::Visit(TotalCrossSection&) const noexcept {}

void Visitor::Visit(TNSL&) const noexcept {}
