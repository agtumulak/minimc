#pragma once

#include "BasicTypes.hpp"
#include "Estimator.hpp"
#include "Particle.hpp"
#include "World.hpp"

#include <vector>

using Bank = std::vector<Particle>;

struct TransportOutcome {
  TransportOutcome& operator+=(const TransportOutcome& rhs) noexcept;
  Estimator estimator;
  Bank banked;
};

/// @brief Try to keep everything a pure function (C++ Core Guidelines F.8)
TransportOutcome TransportAndBank(Bank& source, const World& w);
