#include "ContinuousReaction.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "HDF5DataSet.hpp"
#include "Nuclide.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "Reaction.hpp"
#include "ScalarField.hpp"
#include "ThermalScattering.hpp"
#include "pugixml.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <variant>

// ContinuousReaction

//// public

std::unique_ptr<const ContinuousReaction> ContinuousReaction::Create(
    const pugi::xml_node& reaction_node, const Nuclide& target,
    const std::optional<ThermalScattering>& tnsl) {
  switch (ToReaction(reaction_node.name())) {
  case Reaction::birth:
    assert(false); // this should have been caught by the validator
    return {};
  case Reaction::capture:
    return std::make_unique<const ContinuousCapture>(reaction_node, target);
  case Reaction::scatter:
    return std::make_unique<const ContinuousScatter>(
        reaction_node, target, tnsl);
  case Reaction::fission:
    return std::make_unique<const ContinuousFission>(reaction_node, target);
  case Reaction::leak:
    assert(false); // this should have been caught by the validator
    return {};
  }
}

ContinuousReaction::ContinuousReaction(
    const pugi::xml_node& reaction_node, const Nuclide& target)
    : target{target}, evaluation{reaction_node} {}

ContinuousReaction::~ContinuousReaction() noexcept {}

MicroscopicCrossSection
ContinuousReaction::GetCellMajorant(const Particle& p) const noexcept {
  return evaluation.xs.at(std::get<ContinuousEnergy>(p.GetEnergy()));
}

// ContinuousCapture

ContinuousCapture::ContinuousCapture(
    const pugi::xml_node& capture_node, const Nuclide& target)
    : ContinuousReaction{capture_node, target} {}

void ContinuousCapture::Interact(
    Particle& p, std::vector<Estimator::Proxy>&) const noexcept {
  p.reaction = Reaction::capture;
}

// ContinuousScatter

//// public

ContinuousScatter::ContinuousScatter(
    const pugi::xml_node& scatter_node, const Nuclide& target,
    const std::optional<ThermalScattering>& tnsl)
    : ContinuousReaction{scatter_node.child("xs"), target}, tnsl{tnsl} {}

MicroscopicCrossSection
ContinuousScatter::GetCellMajorant(const Particle& p) const noexcept {
  if (tnsl.has_value() && tnsl->IsValid(p)) {
    return tnsl->GetCellMajorant(p);
  }
  else if (p.type == Particle::Type::neutron) {
    // unadjust evaluated temperature then readjust cross section to requested
    // temperature
    const auto T_max = p.GetCell().temperature->upper_bound;
    const auto free_gas_majorant =
        evaluation.xs.at(std::get<ContinuousEnergy>(p.GetEnergy())) /
        GetFreeGasScatterAdjustment(p, evaluation.temperature) *
        GetFreeGasScatterAdjustment(p, T_max);
    return free_gas_majorant;
  }
  else {
    return ContinuousReaction::GetCellMajorant(p);
  }
}

void ContinuousScatter::Interact(
    Particle& p, std::vector<Estimator::Proxy>&) const noexcept {
  p.reaction = Reaction::scatter;
  if (tnsl.has_value() && tnsl->IsValid(p)) {
    tnsl->Scatter(p);
  }
  else if (p.type == Particle::Type::neutron) {
    // get the true cross section
    const auto T = p.GetCell().temperature->at(p.GetPosition());
    const auto T_max = p.GetCell().temperature->upper_bound;
    // compute probability of a true (non-virtual) collision
    const auto p_true = GetFreeGasScatterAdjustment(p, T) /
                        GetFreeGasScatterAdjustment(p, T_max);
    if (std::bernoulli_distribution{p_true}(p.rng)) {
      // true collision happens; particle state is changed
      ScatterFreeGas(p);
    }
  }
  else {
    std::runtime_error("Scattering for non-neutrons not implemented.");
  }
}

//// private

Real ContinuousScatter::GetFreeGasScatterAdjustment(
    const Particle& p, Temperature T) const noexcept {
  if (T == 0) {
    return 1;
  }
  else {
    const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
    const auto x = std::sqrt(E / (constants::boltzmann * T));
    const auto arg = target.awr * x * x;
    return (1 + 1 / (2 * arg)) * std::erf(std::sqrt(arg)) +
           std::exp(-arg) / std::sqrt(constants::pi * arg);
  }
}

void ContinuousScatter::ScatterFreeGas(Particle& p) const noexcept {
  // Adapted from
  // openmc: openmc.readthedocs.io/en/stable/methods/neutron_physics.html
  // lanl: laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-9721.pdf

  // neutron mass in MeV * (s / cm)^2
  constexpr auto m_n = constants::neutron_mass;
  // incident neutron energy in MeV
  const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
  // neutron speed and velocity in lab frame in cm / s
  const auto s_n = std::sqrt(2. * E / m_n);
  const auto v_n = s_n * p.GetDirection();
  // beta has units of s / cm
  const auto awr = target.awr;
  const auto beta = std::sqrt(
      (awr * m_n) / (2. * constants::boltzmann *
                     p.GetCell().temperature->at(p.GetPosition())));
  // neutron speed (known) and unitless target speed (to be sampled),
  // respectively
  const auto y = beta * s_n;
  Real x;
  // mu is cosine of angle between v_n and target velocity v_T in the lab
  // frame before collision. TODO: Check if this can be made const
  Real mu;
  do {
    // order random numbers are sampled matters: stackoverflow.com/a/40773451
    const auto xi_1 = p.Sample(), xi_2 = p.Sample();
    if (p.Sample() < 2 / (std::sqrt(constants::pi) * y + 2)) { // openmc Eq. 75
      x = std::sqrt(-std::log(xi_1 * xi_2));                   // lanl C49
    }
    else {
      // order random numbers are sampled matters: stackoverflow.com/a/40773451
      const auto xi_3 = p.Sample();
      const auto z = std::cos(constants::pi * xi_3 / 2);
      x = std::sqrt(-std::log(xi_1) - std::log(xi_2) * z * z); // lanl C61
    }
    mu = 2 * p.Sample() - 1; // TODO: Check if this can be factored out
  } while (p.Sample() >= std::sqrt(x * x + y * y - 2 * x * y * mu) / (x + y));
  // the sampled speed of the target in the lab frame
  const auto s_T = x / beta;
  // sample phi, the azimuthal angle about v_n
  const auto phi = 2 * constants::pi * p.Sample();
  // given v_n, s_T, mu, and phi, construct v_T, the velocity of the target
  // in the lab frame.
  const auto v_T = s_T * Direction{p.GetDirection(), mu, phi};
  // Now v_T is known. Compute v_cm, the velocity of the center of mass,
  // followed by V_n, the velocity of the neutron in the CM frame.
  const auto v_cm = (v_n + awr * v_T) / (1 + awr);
  const auto V_n = v_n - v_cm;
  // scattering cosine and azimuthal angle of the neutron in the CM frame
  const auto mu_cm = 2 * p.Sample() - 1;
  const auto phi_cm = 2 * constants::pi * p.Sample();
  // rotate incident neutron velocity to outgoing (primed) neutron velocity
  // in CM
  const auto V_n_prime_direction = Direction{V_n, mu_cm, phi_cm};
  // magnitude of incident velocity V_n is equal to magnitude of outgoing
  // velocity V_n_prime in CM
  const auto V_n_prime = std::sqrt(V_n.Dot(V_n)) * V_n_prime_direction;
  // outgoing velocity in lab frame
  const auto v_n_prime = V_n_prime + v_cm;
  // outgoing energy in lab frame
  const auto E_prime = 0.5 * m_n * v_n_prime.Dot(v_n_prime);
  // update Particle state
  p.SetEnergy(E_prime);
  p.SetDirection(Direction{v_n_prime});
}

// ContinuousFission

//// public

ContinuousFission::ContinuousFission(
    const pugi::xml_node& fission_node, const Nuclide& target)
    : ContinuousReaction{fission_node.child("xs"), target},
      nubar{
          fission_node.child("nubar")
              ? std::make_optional(HDF5DataSet<1>{
                    fission_node.child("nubar").attribute("file").as_string()}
                                       .ToContinuousMap())
              : std::nullopt} {}

void ContinuousFission::Interact(
    Particle& p, std::vector<Estimator::Proxy>&) const noexcept {
  p.reaction = Reaction::fission;
  // rely on the fact that double to int conversions essentially do a floor()
  size_t fission_yield(
      nubar.value().at(std::get<ContinuousEnergy>(p.GetEnergy())) +
      std::uniform_real_distribution{}(p.rng));
  for (size_t i = 0; i < fission_yield; i++) {
    // evaluation order of arguments is undefined so do evaluation here
    const Direction direction{p.rng};
    // todo: add energy sampling
    const auto energy{p.GetEnergy()};
    p.BankSecondaries(direction, energy);
  }
}
