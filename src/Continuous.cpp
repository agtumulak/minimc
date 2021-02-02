#include "Continuous.hpp"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <variant>

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

Continuous::CrossSection
Continuous::GetTotal(const Particle& p) const noexcept {
  return total.at(std::get<ContinuousEnergy>(p.GetEnergy()));
}

void Continuous::Scatter(RNG& rng, Particle& p) const noexcept {
  p.SetDirectionIsotropic(rng);
  return;
}

std::vector<Particle>
Continuous::Fission(RNG& rng, Particle& p) const noexcept {
  std::vector<Particle> fission_neutrons;
  p.Kill();
  // rely on the fact that double to int conversions essentially do a floor()
  size_t fission_yield(
      nubar.value().at(std::get<ContinuousEnergy>(p.GetEnergy())) +
      std::uniform_real_distribution{}(rng));
  for (size_t i = 0; i < fission_yield; i++) {
    // evaluation order of arguments is undefined so do evaluation here
    const auto direction{Direction::CreateIsotropic(rng)};
    const auto energy{Energy{ContinuousEnergy{chi.value().Sample(rng)}}};
    fission_neutrons.emplace_back(
        p.GetPosition(), direction, energy, Particle::Type::neutron);
  }
  assert(fission_neutrons.size() == fission_yield);
  return fission_neutrons;
}

NuclearData::Reaction
Continuous::SampleReaction(RNG& rng, const Particle& p) const noexcept {
  const CrossSection threshold =
      std::uniform_real_distribution{}(rng)*GetTotal(p);
  CrossSection accumulated{0};
  for (const auto& [reaction, xs] : reactions) {
    accumulated += xs.at(std::get<ContinuousEnergy>(p.GetEnergy()));
    if (accumulated > threshold) {
      return reaction;
    }
  }
  // If no reaction found, resample tail-recursively
  return SampleReaction(rng, p);
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
  NuclearData::CrossSection xs;
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
  // Multiply each probability density by bin width. Use previous left value
  // as height as this is a left Riemann sum
  for (auto it = scratch.rbegin(); it != std::prev(scratch.rend()); it++) {
    // multiplication by 1e6 is because energies were stored as MeV and
    // datafile PDF is given in per eV and
    it->second =
        (it->first - std::next(it)->first) * std::next(it)->second * 1e6;
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
    const auto reaction{NuclearData::ToReaction(reaction_name)};
    if (reaction == NuclearData::Reaction::fission) {
      reactions.emplace(
          reaction,
          OneDimensional{
              reaction_node.child("xs").attribute("file").as_string()});
    }
    else {
      reactions.emplace(
          NuclearData::ToReaction(reaction_name),
          OneDimensional{reaction_node.attribute("file").as_string()});
    }
  }
  return reactions;
}
