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
    : nubar{particle_node.child("fission").child("nubar")
              ? std::make_optional(
                  ReadJanisWeb(
                  particle_node.child("fission").child("nubar")
                  .attribute("file").as_string()))
              : std::nullopt},
      chi{particle_node.child("fission").child("chi")
              ? std::make_optional(
                  ReadJanisWebCDF(
                  particle_node.child("fission").child("chi").attribute("file")
                  .as_string()))
              : std::nullopt},
      tsl{ReadPandasSAB(particle_node.child("scatter").child("tsl"))},
      reactions{CreateReactions(particle_node)},
      total{ReadJanisWeb(
          particle_node.child("total").attribute("file").as_string())},
      awr{particle_node.parent().attribute("awr").as_double()}{}

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
           GetAdjustedFreeGasScatter(p, T_max);
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
           GetReaction(p, Reaction::fission) + GetAdjustedFreeGasScatter(p, T);
  }
  else {
    return total.at(std::get<ContinuousEnergy>(p.GetEnergy()));
  }
}

MicroscopicCrossSection
Continuous::GetReaction(const Particle& p, const Reaction r) const noexcept {
  try {
    return reactions.at(r).at(std::get<ContinuousEnergy>(p.GetEnergy()));
  }
  catch (const std::out_of_range& e) {
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

Continuous::CE_XS::elements_type
Continuous::ReadJanisWeb(const std::filesystem::path& datapath) {
  CE_XS::elements_type map;
  std::ifstream datafile{datapath};
  if (!datafile) {
    throw std::runtime_error("File not found: " + datapath.string());
  }
  for (size_t i = 1; i <= 3; i++) {
    // first three lines are header
    datafile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  ContinuousEnergy energy;
  MicroscopicCrossSection xs;
  while (datafile >> energy) {
    datafile.ignore(
        std::numeric_limits<std::streamsize>::max(), ';'); // delimiter
    datafile >> xs;
    map[energy / 1e6] = xs; // datafile is given in eV
  }
  return map;
}

CDF<ContinuousEnergy>::elements_type
Continuous::ReadJanisWebCDF(const std::filesystem::path& datapath) {
  // scratch contains PDF values
  auto scratch = ReadJanisWeb(datapath);
  if (scratch.size() < 2) {
    throw std::runtime_error(
        "In file " + datapath.string() +
        ": at least two entries required to define a CDF");
  }
  // Multiply each probability density by bin width. Use trapezoid rule.
  for (auto it = scratch.rbegin(); it != std::prev(scratch.rend()); it++) {
    // multiplication by 1e6 is because energies were stored as MeV and
    // datafile PDF is given in per eV
    it->second = (it->first - std::next(it)->first) * 0.5 *
                 (it->second + std::next(it)->second) * 1e6;
  }
  // set first element to zero to make accumulation pass look clean af
  scratch.begin()->second = 0;
  for (auto it = std::next(scratch.begin()); it != scratch.end(); it++) {
    it->second += std::prev(it)->second;
  }
  scratch.erase(scratch.cbegin());
  // normalize CDF
  const auto total_weight = scratch.crbegin()->second;
  if (total_weight == 0.) {
    throw std::runtime_error(
        "In file " + datapath.string() + ": total weight is zero.");
  }
  std::for_each(scratch.begin(), scratch.end(), [&total_weight](auto& element) {
    element.second /= total_weight;
  });
  // swap keys and values, this will remove any zero probability values
  CDF<ContinuousEnergy>::elements_type swapped;
  for (const auto& element : scratch) {
    swapped.emplace(element.second, element.first);
  }
  return swapped;
}

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

std::map<Reaction, Continuous::CE_XS>
Continuous::CreateReactions(const pugi::xml_node& particle_node) {
  std::map<Reaction, Continuous::CE_XS> reactions;
  for (const auto& reaction_node : particle_node) {
    const std::string reaction_name = reaction_node.name();
    if (reaction_name == "total") {
      continue; // skip total cross section
    }
    const auto reaction{ToReaction(reaction_name)};
    if (reaction == Reaction::capture) {
      reactions.emplace(
          ToReaction(reaction_name),
          ReadJanisWeb(reaction_node.attribute("file").as_string()));
    }
    else {
      reactions.emplace(
          reaction,
          ReadJanisWeb(
              reaction_node.child("xs").attribute("file").as_string()));
    }
  }
  return reactions;
}

void Continuous::Capture(Particle& p) const noexcept { p.Kill(); }

void Continuous::Scatter(Particle& p) const noexcept {
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
      if (p.Sample() <
          2 / (std::sqrt(constants::pi) * y + 2)) {        // openmc Eq. 75
        x = std::sqrt(-std::log(p.Sample() * p.Sample())); // lanl C49
      }
      else {
        // order random numbers are sampled matters: stackoverflow.com/a/40773451
        const auto xi_1 = p.Sample(), xi_2 = p.Sample(), xi_3 = p.Sample();
        const auto z = std::cos(constants::pi * xi_3 / 2);
        x = std::sqrt(-std::log(xi_1) - std::log(xi_2) * z * z); // lanl C61
      }
      mu = 2 * p.Sample() - 1; // TODO: Check if this can be factored out
    } while (p.Sample() < std::sqrt(x * x + y * y - 2 * x * y * mu) / (x + y));
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
    // Compute neutron energy in CM frame. Use classical kinetic energy. Speed
    // units are cm / s. Energy units are MeV.
    const auto E_n = 0.5 * m_n * V_n.Dot(V_n);
    // scattering cosine of the neutron in the CM frame
    const auto mu_cm = 2 * p.Sample() - 1;
    // Compute e_n, the outgoing energy in lab frame. From openmc Eq. 52
    const auto e_n = E_n + (E + 2 * mu_cm * (awr + 1) * std::sqrt(E * E_n)) /
                               ((awr + 1) * (awr + 1));
    // mu_lab is the scattering cosine in the lab frame.
    const auto mu_lab =
        mu * std::sqrt(E_n / e_n) + 1 / (awr + 1) * std::sqrt(E / e_n);
    // Update particle state
    p.Scatter(mu_lab, e_n);
  }
  return;
}

void Continuous::Fission(Particle& p) const noexcept {
  p.Kill();
  // rely on the fact that double to int conversions essentially do a floor()
  size_t fission_yield(
      nubar.value().at(std::get<ContinuousEnergy>(p.GetEnergy())) +
      std::uniform_real_distribution{}(p.rng));
  for (size_t i = 0; i < fission_yield; i++) {
    // evaluation order of arguments is undefined so do evaluation here
    const auto direction{Direction::CreateIsotropic(p.rng)};
    const auto energy{Energy{ContinuousEnergy{chi.value().Sample(p.rng)}}};
    p.Bank(direction, energy);
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

MicroscopicCrossSection Continuous::GetAdjustedFreeGasScatter(
    const Particle& p, Temperature T) const noexcept {
  const auto E = std::get<ContinuousEnergy>(p.GetEnergy());
  const auto x = std::sqrt(E / (constants::boltzmann * T));
  const auto arg = awr * x * x;
  return GetReaction(p, Reaction::scatter) *
         ((1 + 1 / (2 * arg)) * std::erf(std::sqrt(arg)) +
          std::exp(-arg) / std::sqrt(constants::pi * arg));
}
