#include "Continuous.hpp"

#include "HDF5DataSet.hpp"
#include "Particle.hpp"
#include "Point.hpp"

#include <algorithm>
#include <cstddef>
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

ContinuousMap<ContinuousEnergy, MicroscopicCrossSection>
Continuous::ReadPandasHDF5(const std::filesystem::path& datapath) {
  std::map<ContinuousEnergy, MicroscopicCrossSection> map;
  HDF5DataSet h5_data{datapath};
  const auto& energies = h5_data.GetAxis(0);
  for (size_t index = 0; index < energies.size(); index++) {
    map[energies[index]] = h5_data.at(index);
  }
  return map;
}

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
          particle_node.child("total").attribute("file").as_string())} {}

MicroscopicCrossSection Continuous::GetTotal(const Particle& p) const noexcept {
  if (tsl->IsValid(p)) {
    return GetReaction(p, Reaction::capture) + tsl->GetMajorant(p);
  }
  return total.at(std::get<ContinuousEnergy>(p.GetEnergy()));
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
  if (tsl->IsValid(p)) {
    tsl->Scatter(p);
  }
  else {
    p.SetDirectionIsotropic();
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
