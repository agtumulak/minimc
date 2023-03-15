#include "Perturbation/Sensitivity/Proxy/Visitor.hpp"

using namespace Perturbation::Sensitivity::Proxy;

// Visitor

//// public

Visitor::Visitor(const Particle& p, const BinIndex i, const Score s) noexcept
    : particle{p}, index{i}, score{s} {}

Visitor::~Visitor() noexcept {}
