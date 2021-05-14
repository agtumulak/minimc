#include "Continuous.hpp"

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
#include <vector>

// Continuous

//// public

Continuous::Continuous(const pugi::xml_node& particle_node)
    : nubar{particle_node.child("fission").child("nubar")
              ? std::make_optional<OneDimensional>(
                  particle_node.child("fission").child("nubar")
                  .attribute("file").as_string())
              : std::nullopt},
      chi{particle_node.child("fission").child("chi")
              ? std::make_optional<CDF>(
                  particle_node.child("fission").child("chi").attribute("file")
                  .as_string())
              : std::nullopt},
      reactions{CreateReactions(particle_node)},
      total{particle_node.child("total").attribute("file").as_string()} {}

MicroscopicCrossSection Continuous::GetTotal(const Particle& p) const noexcept {
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

// Continuous::OneDimensional

//// public

Continuous::OneDimensional::OneDimensional(
    const std::filesystem::path& datapath) {
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
    elements[energy / 1e6] = xs; // datafile is given in eV
  }
}

const Real&
Continuous::OneDimensional::at(const ContinuousEnergy e) const noexcept {
  // TODO: Interpolate and handle edge cases
  return (*elements.upper_bound(e)).second;
}

// Continuous::CDF

//// public

Continuous::CDF::CDF(const std::filesystem::path& datapath)
    : OneDimensional{datapath} {
  elements_type scratch = std::move(elements);
  elements.clear(); // the final map will use CDF values as keys
  if (scratch.size() < 2) {
    throw std::runtime_error(
        "In file " + datapath.string() +
        ": at least two entries required to define a CDF");
  }
  // Multiply each probability density by bin width. Use trapezoid rule.
  for (auto it = scratch.rbegin(); it != std::prev(scratch.rend()); it++) {
    // multiplication by 1e6 is because energies were stored as MeV and
    // datafile PDF is given in per eV and
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
  for (const auto& element : scratch) {
    elements.emplace(element.second, element.first);
  }
}

Continuous::elements_type::mapped_type
Continuous::CDF::Sample(RNG& rng) const noexcept {
  return elements.upper_bound(std::uniform_real_distribution{}(rng))->second;
}

//  Continuous

//// private

Continuous::ReactionsMap
Continuous::CreateReactions(const pugi::xml_node& particle_node) {
  Continuous::ReactionsMap reactions;
  for (const auto& reaction_node : particle_node) {
    const std::string reaction_name = reaction_node.name();
    if (reaction_name == "total") {
      continue; // skip total cross section
    }
    const auto reaction{ToReaction(reaction_name)};
    if (reaction == Reaction::fission) {
      reactions.emplace(
          reaction,
          OneDimensional{
              reaction_node.child("xs").attribute("file").as_string()});
    }
    else {
      reactions.emplace(
          ToReaction(reaction_name),
          OneDimensional{reaction_node.attribute("file").as_string()});
    }
  }
  return reactions;
}

void Continuous::Capture(Particle& p) const noexcept { p.Kill(); }

void Continuous::Scatter(Particle& p) const noexcept {
  p.SetDirectionIsotropic();
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
    p.secondaries.emplace_back(
        p.GetPosition(), direction, energy, Particle::Type::neutron);
    p.secondaries.back().SetCell(p.GetCell());
    p.secondaries.back().seed =
        std::uniform_int_distribution<RNG::result_type>{1}(p.rng);
  }
}
