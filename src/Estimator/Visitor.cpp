#include "Estimator/Visitor.hpp"

using namespace Estimator;

// EstimatorVisitor

//// public

Visitor::Visitor(const Particle& p) noexcept : particle{p} {};

Visitor::~Visitor() noexcept {}
