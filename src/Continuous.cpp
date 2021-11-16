#include "Continuous.hpp"

#include "Cell.hpp"
#include "ContinuousMap.hpp"
#include "ContinuousReaction.hpp"
#include "Particle.hpp"
#include "ScalarField.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <iosfwd>
#include <numeric>
#include <random>
#include <string>
#include <type_traits>
#include <variant>

// Continuous

//// public

Continuous::Continuous(const pugi::xml_node& particle_node)
    : reactions{CreateReactions(particle_node)}, total{particle_node.child(
                                                     "total")} {}

MicroscopicCrossSection
Continuous::GetMajorant(const State& s) const noexcept {
  // majorant temperature is assumed to occur at maximum temperature in Cell
  const auto majorant_temperature = s.cell->temperature->upper_bound;
  if (!ReactionsModifyTotal(s) && total.IsValid(majorant_temperature)) {
    return total.xs.at(std::get<ContinuousEnergy>(s.energy));
  }
  else {
    return std::accumulate(
        reactions.cbegin(), reactions.cend(), MicroscopicCrossSection{0},
        [&s](const auto& accumulated, const auto& reaction) {
          return accumulated + reaction->GetMajorant(s);
        });
  }
}

MicroscopicCrossSection Continuous::GetTotal(const State& s) const noexcept {
  const auto requested_temperature = s.cell->temperature->at(s.position);
  if (!ReactionsModifyTotal(s) && total.IsValid(requested_temperature)) {
    return total.xs.at(std::get<ContinuousEnergy>(s.energy));
  }
  else {
    return std::accumulate(
        reactions.cbegin(), reactions.cend(), MicroscopicCrossSection{0},
        [&s](const auto& accumulated, const auto& reaction) {
          return accumulated + reaction->GetCrossSection(s);
        });
  }
}

void Continuous::Interact(State& s) const noexcept {
  const MicroscopicCrossSection threshold =
      std::uniform_real_distribution{}(s.rng) * GetTotal(s);
  MicroscopicCrossSection accumulated{0};
  for (const auto& candidate : reactions) {
    accumulated += candidate->GetCrossSection(s);
    if (accumulated > threshold) {
      candidate->Interact(s);
      return;
    }
  }
  // If no reaction found, resample tail-recursively
  return Interact(s);
}

//// private

std::vector<std::unique_ptr<const ContinuousReaction>>
Continuous::CreateReactions(const pugi::xml_node& particle_node) {
  std::vector<std::unique_ptr<const ContinuousReaction>> reactions;
  for (const auto& reaction_node : particle_node) {
    const std::string reaction_name = reaction_node.name();
    if (reaction_name == "total") {
      continue; // skip total cross section
    }
    reactions.emplace_back(ContinuousReaction::Create(reaction_node));
  }
  return reactions;
}

bool Continuous::ReactionsModifyTotal(const State& s) const noexcept {
  return std::any_of(
      reactions.cbegin(), reactions.cend(),
      [&s](const auto& reaction) { return reaction->ModifiesTotal(s); });
}
