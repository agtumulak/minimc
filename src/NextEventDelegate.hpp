#pragma once

#include "Particle.hpp"

class NextEventDelegate {
public:
  virtual void Sample(Particle& p) const noexcept = 0;
};

class SurfaceTracking : public NextEventDelegate {
public:
  void Sample(Particle& p) const noexcept override;
};
