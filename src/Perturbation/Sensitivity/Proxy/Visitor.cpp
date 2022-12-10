#include "Perturbation/Sensitivity/Proxy/Visitor.hpp"

using namespace Perturbation::Sensitivity::Proxy;

// Visitor

//// public

Visitor::Visitor(
    const Particle& p, const BinIndex i, const Score s,
    const bool add_to_existing) noexcept
    : particle{p}, index{i}, score{s}, add_to_existing{add_to_existing} {}

Visitor::~Visitor() noexcept {}
