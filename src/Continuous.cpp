#include "Continuous.hpp"

#include "Cell.hpp"
#include "Constants.hpp"
#include "HDF5DataSet.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "ScalarField.hpp"

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>

// Continuous

//// public

Continuous::Continuous(const pugi::xml_node& particle_node)
    : nubar{
      particle_node.child("fission").child("nubar")
        ? std::make_optional(
            HDF5DataSet<1>{
            particle_node.child("fission").child("nubar").attribute("file")
            .as_string()}.ToContinuousMap())
        : std::nullopt},
      tsl{ReadPandasSAB(particle_node.child("scatter").child("tsl"))},
      reactions{CreateReactions(particle_node)},
      total{HDF5DataSet<1>{
          particle_node.child("total").attribute("file").as_string()}
                .ToContinuousMap()},
      awr{particle_node.parent().attribute("awr").as_double()} {}

MicroscopicCrossSection
Continuous::GetMajorant(const Particle& p) const noexcept {
  if (tsl.has_value() && tsl->IsValid(p)) {
    // thermal scattering is assumed to encompass elastic and inelstic
    // scattering
    // TODO: Possibly make this polymorphic to avoid `if`?
    return GetReaction(p, Reaction::capture) +
           GetReaction(p, Reaction::fission) + tsl->GetMajorant(p);
  }
  else if (const auto T_max = p.GetCell().temperature->upper_bound;
           IsFreeGasScatteringValid(p, T_max)) {
    return GetReaction(p, Reaction::capture) +
           GetReaction(p, Reaction::fission) +
           GetReaction(p, Reaction::scatter) *
               GetFreeGasScatterAdjustment(p, T_max);
  }
  else {
    return GetTotal(p);
  }
}

MicroscopicCrossSection Continuous::GetTotal(const Particle& p) const noexcept {
  if (tsl.has_value() && tsl->IsValid(p)){
    return GetReaction(p, Reaction::capture) +
           GetReaction(p, Reaction::fission) + tsl->GetTotal(p);
  }
  else if (const auto T = p.GetCell().temperature->at(p.GetPosition());
           IsFreeGasScatteringValid(p, T)) {
    return GetReaction(p, Reaction::capture) +
           GetReaction(p, Reaction::fission) +
           GetReaction(p, Reaction::scatter) *
               GetFreeGasScatterAdjustment(p, T);
  }
  else {
    return total.at(std::get<ContinuousEnergy>(p.GetEnergy()));
  }
}

MicroscopicCrossSection
Continuous::GetReaction(const Particle& p, const Reaction r) const noexcept {
  if (const auto reaction_it = reactions.find(r);
      reaction_it != reactions.cend()) {
    return reaction_it->second.at(std::get<ContinuousEnergy>(p.GetEnergy()));
  }
  else {
    return 0;
  }
}

Real Continuous::GetNuBar(const Particle& p) const noexcept {
  if (nubar) {
    return nubar->at(std::get<ContinuousEnergy>(p.GetEnergy()));
  }
  else {
    return 0;
  }
}

void Continuous::Interact(Particle& p) const noexcept {
  const MicroscopicCrossSection threshold =
      std::uniform_real_distribution{}(p.rng) * GetTotal(p);
  MicroscopicCrossSection accumulated{0};
  for (const auto& [reaction, xs] : reactions) {
    accumulated += xs.at(std::get<ContinuousEnergy>(p.GetEnergy()));
    if (accumulated > threshold) {
      switch (reaction) {
      case Reaction::capture:
        Capture(p);
        break;
      case Reaction::scatter:
        Scatter(p);
        break;
      case Reaction::fission:
        Fission(p);
        break;
      }
      return;
    }
  }
  // If no reaction found, resample tail-recursively
  return Interact(p);
}

//// private

std::optional<ThermalScattering>
Continuous::ReadPandasSAB(const pugi::xml_node& tsl_node) {
  if (!tsl_node) {
    return std::nullopt;
  }
  if (Particle::ToType(tsl_node.parent().parent().name()) !=
      Particle::Type::neutron) {
    throw std::runtime_error(
        tsl_node.path() +
        ": Only neutrons may have a thermal scattering library node");
  }
  return ThermalScattering{tsl_node};
}

std::map<Reaction, ContinuousMap<ContinuousEnergy, MicroscopicCrossSection>>
Continuous::CreateReactions(const pugi::xml_node& particle_node) {
  std::map<Reaction, ContinuousMap<ContinuousEnergy, MicroscopicCrossSection>>
      reactions;
  for (const auto& reaction_node : particle_node) {
    const std::string reaction_name = reaction_node.name();
    if (reaction_name == "total") {
      continue; // skip total cross section
    }
    const auto reaction{ToReaction(reaction_name)};
    if (reaction == Reaction::capture) {
      reactions.emplace(
          reaction, HDF5DataSet<1>{reaction_node.attribute("file").as_string()}
                        .ToContinuousMap());
    }
    else {
      reactions.emplace(
          reaction,
          HDF5DataSet<1>{
              reaction_node.child("xs").attribute("file").as_string()}
              .ToContinuousMap());
    }
  }
  return reactions;
}

void Continuous::Capture(Particle& p) const noexcept {
  p.event = Particle::Event::capture;
}

void Continuous::Scatter(Particle& p) const noexcept {
  p.event = Particle::Event::scatter;
  if (tsl.has_value() && tsl->IsValid(p)) {
    tsl->Scatter(p);
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
          2 / (std::sqrt(constants::pi) * y + 2)) {        // openmc Eq. 75
        x = std::sqrt(-std::log(xi_1 * xi_2)); // lanl C49
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
    const auto v_T =
        s_T * Direction::CreateAboutDirection(p.GetDirection(), mu, phi);
    // Now v_T is known. Compute v_cm, the velocity of the center of mass,
    // followed by V_n, the velocity of the neutron in the CM frame.
    const auto v_cm = (v_n + awr * v_T) / (1 + awr);
    const auto V_n = v_n - v_cm;
    // scattering cosine and azimuthal angle of the neutron in the CM frame
    const auto mu_cm = 2 * p.Sample() - 1;
    const auto phi_cm = 2 * constants::pi * p.Sample();
    // rotate incident neutron velocity to outgoing (primed) neutron velocity
    // in CM
    const auto V_n_prime_direction =
        Direction::CreateAboutDirection(V_n, mu_cm, phi_cm);
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

void Continuous::Fission(Particle& p) const noexcept {
  p.event = Particle::Event::fission;
  // rely on the fact that double to int conversions essentially do a floor()
  size_t fission_yield(
      nubar.value().at(std::get<ContinuousEnergy>(p.GetEnergy())) +
      std::uniform_real_distribution{}(p.rng));
  for (size_t i = 0; i < fission_yield; i++) {
    // evaluation order of arguments is undefined so do evaluation here
    const auto direction{Direction::CreateIsotropic(p.rng)};
    // todo: add energy sampling
    const auto energy{p.GetEnergy()};
    p.BankSecondaries(direction, energy);
  }
}

bool Continuous::IsFreeGasScatteringValid(
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

Real Continuous::GetFreeGasScatterAdjustment(
    const Particle& p, Temperature T) const noexcept {
  const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
  const auto x = std::sqrt(E / (constants::boltzmann * T));
  const auto arg = awr * x * x;
  return (1 + 1 / (2 * arg)) * std::erf(std::sqrt(arg)) +
         std::exp(-arg) / std::sqrt(constants::pi * arg);
}
