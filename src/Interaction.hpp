#pragma once

#include "BasicTypes.hpp"

class Particle;

/// @brief Models the interaction between a Particle and a Nuclide
/// @details The polymorphism here shall be where multigroup and continuous
///          energy cross sections are resolved.
class Interaction {
public:
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Interaction() noexcept;
  /// @brief Returns the majorant cross section for a given Particle
  virtual MicroscopicCrossSection
  GetMajorant(const Particle& p) const noexcept = 0;
  /// @brief Returns the total cross section for a given Particle
  virtual MicroscopicCrossSection
  GetTotal(const Particle& p) const noexcept = 0;
  /// @brief Interact with a Particle, updating its state
  virtual void Interact(Particle& p) const noexcept = 0;
};
