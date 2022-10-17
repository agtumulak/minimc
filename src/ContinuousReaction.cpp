#include "ContinuousReaction.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "HDF5DataSet.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "Reaction.hpp"
#include "ScalarField.hpp"
#include "pugixml.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>
#include <iostream>

// ContinuousReaction

//// public

std::unique_ptr<const ContinuousReaction>
ContinuousReaction::Create(const pugi::xml_node& reaction_node) {
  switch (ToReaction(reaction_node.name())) {
  case Reaction::capture:
    return std::make_unique<const ContinuousCapture>(reaction_node);
  case Reaction::scatter:
    return std::make_unique<const ContinuousScatter>(reaction_node);
  case Reaction::fission:
    return std::make_unique<const ContinuousFission>(reaction_node);
  }
}

ContinuousReaction::ContinuousReaction(const pugi::xml_node& reaction_node)
    : evaluation{reaction_node} {}

ContinuousReaction::~ContinuousReaction() noexcept {}

bool ContinuousReaction::ModifiesTotal(const Particle&) const noexcept {
  return false;
}

MicroscopicCrossSection
ContinuousReaction::GetMajorant(const Particle& p) const noexcept {
  return GetCrossSection(p);
}

MicroscopicCrossSection
ContinuousReaction::GetCrossSection(const Particle& p) const noexcept {
  return evaluation.xs.at(std::get<ContinuousEnergy>(p.GetEnergy()));
}

// ContinuousCapture

ContinuousCapture::ContinuousCapture(const pugi::xml_node& capture_node)
    : ContinuousReaction{capture_node} {}

void ContinuousCapture::Interact(Particle& p) const noexcept {
  p.event = Particle::Event::capture;
}

// ContinuousScatter

//// public

ContinuousScatter::ContinuousScatter(const pugi::xml_node& scatter_node)
    : ContinuousReaction{scatter_node.child("xs")},
      tnsl{ReadPandasSAB(scatter_node.child("tnsl"))},
      awr{scatter_node.parent().parent().attribute("awr").as_double()} {}

bool ContinuousScatter::ModifiesTotal(const Particle& p) const noexcept {
  return tnsl.has_value() && tnsl->IsValid(p);
}

MicroscopicCrossSection
ContinuousScatter::GetMajorant(const Particle& p) const noexcept {
  if (tnsl.has_value() && tnsl->IsValid(p)) {
    return tnsl->GetMajorant(p);
  }
  else if (const auto T_max = p.GetCell().temperature->upper_bound;
           evaluation.IsValid(T_max)) {
    return evaluation.xs.at(std::get<ContinuousEnergy>(p.GetEnergy()));
  }
  else if (IsFreeGasScatteringValid(p, T_max)) {
    // unadjust evaluated temperature then readjust cross section to requested
    // temperature
    return evaluation.xs.at(std::get<ContinuousEnergy>(p.GetEnergy())) /
           GetFreeGasScatterAdjustment(p, evaluation.temperature) *
           GetFreeGasScatterAdjustment(p, T_max);
  }
  else {
    return evaluation.xs.at(std::get<ContinuousEnergy>(p.GetEnergy()));
  }
}

MicroscopicCrossSection
ContinuousScatter::GetCrossSection(const Particle& p) const noexcept {
  if (tnsl.has_value() && tnsl->IsValid(p)) {
    return tnsl->GetTotal(p);
  }
  else if (const auto T = p.GetCell().temperature->at(p.GetPosition());
           evaluation.IsValid(T)) {
    return evaluation.xs.at(std::get<ContinuousEnergy>(p.GetEnergy()));
  }
  else if (IsFreeGasScatteringValid(p, T)) {
    // unadjust evaluated temperature then readjust cross section to requested
    // temperature
    return evaluation.xs.at(std::get<ContinuousEnergy>(p.GetEnergy())) /
           GetFreeGasScatterAdjustment(p, evaluation.temperature) *
           GetFreeGasScatterAdjustment(p, T);
  }
  else {
    return evaluation.xs.at(std::get<ContinuousEnergy>(p.GetEnergy()));
  }
}

void ContinuousScatter::Interact(Particle& p) const noexcept {
  p.event = Particle::Event::scatter;
  if (tnsl.has_value() && tnsl->IsValid(p)) {
    tnsl->Scatter(p);
  }
  else {
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
      if (p.Sample() <
          2 / (std::sqrt(constants::pi) * y + 2)) { // openmc Eq. 75
        x = std::sqrt(-std::log(xi_1 * xi_2));      // lanl C49
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
    p.Scatter(p.GetDirection().Dot(Direction{v_n_prime}), E_prime);
  }
  return;
}

//// private

std::optional<ThermalScattering>
ContinuousScatter::ReadPandasSAB(const pugi::xml_node& tnsl_node) {
  if (!tnsl_node) {
    return std::nullopt;
  }
  if (Particle::ToType(tnsl_node.parent().parent().name()) !=
      Particle::Type::neutron) {
    throw std::runtime_error(
        tnsl_node.path() +
        ": Only neutrons may have a thermal scattering library node");
  }
  return ThermalScattering{tnsl_node};
}

bool ContinuousScatter::IsFreeGasScatteringValid(
    const Particle& p, const Temperature& T) const noexcept {
  if (p.type != Particle::Type::neutron) {
    return false;
  }
  else if (awr <= 1.0) {
    return true; // hydrogen must always be treated with free gas scattering
  }
  else if (
      std::get<ContinuousEnergy>(p.GetEnergy()) <
      500 * constants::boltzmann * T / awr) {
    return true;
  }
  else {
    return false;
  }
}

Real ContinuousScatter::GetFreeGasScatterAdjustment(
    const Particle& p, Temperature T) const noexcept {
  if (T == 0) {
    return 1;
  }
  else {
    const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
    const auto x = std::sqrt(E / (constants::boltzmann * T));
    const auto arg = awr * x * x;
    return (1 + 1 / (2 * arg)) * std::erf(std::sqrt(arg)) +
           std::exp(-arg) / std::sqrt(constants::pi * arg);
  }
}

// ContinuousFission

//// public

ContinuousFission::ContinuousFission(const pugi::xml_node& fission_node)
    : ContinuousReaction{fission_node.child("xs")},
      nubar{
          fission_node.child("nubar")
              ? std::make_optional(HDF5DataSet<1>{
                    fission_node.child("nubar").attribute("file").as_string()}
                                       .ToContinuousMap())
              : std::nullopt} {}

void ContinuousFission::Interact(Particle& p) const noexcept {
  p.event = Particle::Event::fission;
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
