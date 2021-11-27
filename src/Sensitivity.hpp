#pragma once

#include "BasicTypes.hpp"
#include "Perturbation.hpp"

namespace pugi {
class xml_node;
}
class CurrentEstimator;
class Estimator;
class Particle;


/// @brief Models the change in an Estimator with respect to a Perturbation
class Sensitivity {
public:
  /// @brief Factory method to create a new Sensitivity from a sensitivity node
  ///        of an XML document
  static std::unique_ptr<Sensitivity>
  Create(const pugi::xml_node& sensitivity_node) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~Sensitivity() noexcept;

  virtual Real GetDirectEffect(const Particle& p) const noexcept = 0;

  virtual Real GetMultiplier(const Particle& p) const noexcept = 0;

private:
  // The Estimator which will respond to a perturbation
  const Estimator& e;
};

/// @brief Used for unmultiplied (the actual score) scoring
class NoSensitivity : public Sensitivity {
public:
  Real GetDirectEffect(const Particle&) const noexcept override;
};

class TotalCrossSectionPerturbation : public Sensitivity {
public:
  /// @brief Returns the indirect and direct effect on a given Estimator
  Real GetDirectEffect(const Particle& p) const noexcept override;
};

class CurrentWrtTotal : public Sensitivity {
public:
  Real GetDirectEffect(const Particle& p) const noexcept override;

private:
  const CurrentEstimator& current_estimator;
};
