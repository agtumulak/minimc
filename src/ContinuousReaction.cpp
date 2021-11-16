#include "ContinuousReaction.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "HDF5DataSet.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "Reaction.hpp"
#include "ScalarField.hpp"
#include "State.hpp"
#include "pugixml.hpp"

#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <variant>

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

bool ContinuousReaction::ModifiesTotal(const State&) const noexcept {
  return false;
}

MicroscopicCrossSection
ContinuousReaction::GetMajorant(const State& s) const noexcept {
  return GetCrossSection(s);
}

MicroscopicCrossSection
ContinuousReaction::GetCrossSection(const State& s) const noexcept {
  return evaluation.xs.at(std::get<ContinuousEnergy>(s.energy));
}

// ContinuousCapture

ContinuousCapture::ContinuousCapture(const pugi::xml_node& capture_node)
    : ContinuousReaction{capture_node} {}

void ContinuousCapture::Interact(State& s) const noexcept {
  s.event = State::Event::capture;
}

// ContinuousScatter

//// public

ContinuousScatter::ContinuousScatter(const pugi::xml_node& scatter_node)
    : ContinuousReaction{scatter_node.child("xs")},
      tsl{ReadPandasSAB(scatter_node.child("tsl"))},
      awr{scatter_node.parent().parent().attribute("awr").as_double()} {}

bool ContinuousScatter::ModifiesTotal(const State& s) const noexcept {
  return tsl.has_value() && tsl->IsValid(s);
}

MicroscopicCrossSection
ContinuousScatter::GetMajorant(const State& s) const noexcept {
  if (tsl.has_value() && tsl->IsValid(s)) {
    return tsl->GetMajorant(s);
  }
  else if (const auto T_max = s.cell->temperature->upper_bound;
           evaluation.IsValid(T_max)) {
    return evaluation.xs.at(std::get<ContinuousEnergy>(s.energy));
  }
  else if (IsFreeGasScatteringValid(s, T_max)) {
    // unadjust evaluated temperature then readjust cross section to requested
    // temperature
    return evaluation.xs.at(std::get<ContinuousEnergy>(s.energy)) /
           GetFreeGasScatterAdjustment(s, evaluation.temperature) *
           GetFreeGasScatterAdjustment(s, T_max);
  }
  else {
    return evaluation.xs.at(std::get<ContinuousEnergy>(s.energy));
  }
}

MicroscopicCrossSection
ContinuousScatter::GetCrossSection(const State& s) const noexcept {
  if (tsl.has_value() && tsl->IsValid(s)) {
    return tsl->GetTotal(s);
  }
  else if (const auto T = s.cell->temperature->at(s.position);
           evaluation.IsValid(T)) {
    return evaluation.xs.at(std::get<ContinuousEnergy>(s.energy));
  }
  else if (IsFreeGasScatteringValid(s, T)) {
    // unadjust evaluated temperature then readjust cross section to requested
    // temperature
    return evaluation.xs.at(std::get<ContinuousEnergy>(s.energy)) /
           GetFreeGasScatterAdjustment(s, evaluation.temperature) *
           GetFreeGasScatterAdjustment(s, T);
  }
  else {
    return evaluation.xs.at(std::get<ContinuousEnergy>(s.energy));
  }
}

void ContinuousScatter::Interact(State& s) const noexcept {
  s.event = State::Event::scatter;
  if (tsl.has_value() && tsl->IsValid(s)) {
    tsl->Scatter(s);
  }
  else {
    // Adapted from
    // openmc: openmc.readthedocs.io/en/stable/methods/neutron_physics.html
    // lanl: laws.lanl.gov/vhosts/mcnp.lanl.gov/pdf_files/la-9721.pdf

    // neutron mass in MeV * (s / cm)^2
    constexpr auto m_n = constants::neutron_mass;
    // incident neutron energy in MeV
    const auto E = std::get<ContinuousEnergy>(s.energy);
    // neutron speed and velocity in lab frame in cm / s
    const auto s_n = std::sqrt(2. * E / m_n);
    const auto v_n = s_n * s.direction;
    // beta has units of s / cm
    const auto beta = std::sqrt(
        (awr * m_n) /
        (2. * constants::boltzmann * s.cell->temperature->at(s.position)));
    // neutron speed (known) and unitless target speed (to be sampled),
    // respectively
    const auto y = beta * s_n;
    Real x;
    // mu is cosine of angle between v_n and target velocity v_T in the lab
    // frame before collision. TODO: Check if this can be made const
    Real mu;
    do {
      // order random numbers are sampled matters: stackoverflow.com/a/40773451
      const auto xi_1 = s.Sample(), xi_2 = s.Sample();
      if (s.Sample() <
          2 / (std::sqrt(constants::pi) * y + 2)) { // openmc Eq. 75
        x = std::sqrt(-std::log(xi_1 * xi_2));      // lanl C49
      }
      else {
        // order random numbers are sampled matters: stackoverflow.com/a/40773451
        const auto xi_3 = s.Sample();
        const auto z = std::cos(constants::pi * xi_3 / 2);
        x = std::sqrt(-std::log(xi_1) - std::log(xi_2) * z * z); // lanl C61
      }
      mu = 2 * s.Sample() - 1; // TODO: Check if this can be factored out
    } while (s.Sample() >= std::sqrt(x * x + y * y - 2 * x * y * mu) / (x + y));
    // the sampled speed of the target in the lab frame
    const auto s_T = x / beta;
    // sample phi, the azimuthal angle about v_n
    const auto phi = 2 * constants::pi * s.Sample();
    // given v_n, s_T, mu, and phi, construct v_T, the velocity of the target
    // in the lab frame.
    const auto v_T = s_T * Direction(s.direction, mu, phi);
    // Now v_T is known. Compute v_cm, the velocity of the center of mass,
    // followed by V_n, the velocity of the neutron in the CM frame.
    const auto v_cm = (v_n + awr * v_T) / (1 + awr);
    const auto V_n = v_n - v_cm;
    // scattering cosine and azimuthal angle of the neutron in the CM frame
    const auto mu_cm = 2 * s.Sample() - 1;
    const auto phi_cm = 2 * constants::pi * s.Sample();
    // rotate incident neutron velocity to outgoing (primed) neutron velocity
    // in CM
    const auto V_n_prime_direction = Direction(V_n, mu_cm, phi_cm);
    // magnitude of incident velocity V_n is equal to magnitude of outgoing
    // velocity V_n_prime in CM
    const auto V_n_prime = std::sqrt(V_n.Dot(V_n)) * V_n_prime_direction;
    // outgoing velocity in lab frame
    const auto v_n_prime = V_n_prime + v_cm;
    // outgoing energy in lab frame
    const auto E_prime = 0.5 * m_n * v_n_prime.Dot(v_n_prime);
    // modify State
    s.direction = Direction{v_n_prime};
    s.energy = E_prime;
  }
  return;
}

//// private

std::optional<ThermalScattering>
ContinuousScatter::ReadPandasSAB(const pugi::xml_node& tsl_node) {
  if (!tsl_node) {
    return std::nullopt;
  }
  if (ToParticle(tsl_node.parent().parent().name()) != Particle::neutron) {
    throw std::runtime_error(
        tsl_node.path() +
        ": Only neutrons may have a thermal scattering library node");
  }
  return ThermalScattering{tsl_node};
}

bool ContinuousScatter::IsFreeGasScatteringValid(
    const State& s, const Temperature& T) const noexcept {
  if (s.particle != Particle::neutron) {
    return false;
  }
  else if (awr <= 1.0) {
    return true; // hydrogen must always be treated with free gas scattering
  }
  else if (
      std::get<ContinuousEnergy>(s.energy) <
      500 * constants::boltzmann * T / awr) {
    return true;
  }
  else {
    return false;
  }
}

Real ContinuousScatter::GetFreeGasScatterAdjustment(
    const State& s, Temperature T) const noexcept {
  if (T == 0) {
    return 1;
  }
  else {
    const auto E = std::get<ContinuousEnergy>(s.energy);
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

void ContinuousFission::Interact(State& s) const noexcept {
  s.event = State::Event::fission;
}
