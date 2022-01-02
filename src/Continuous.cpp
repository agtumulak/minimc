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
Continuous::GetMajorant(const Particle& p) const noexcept {
  // majorant cross section is assumed to occur at maximum temperature in Cell
  const auto majorant_temperature = p.GetCell().temperature->upper_bound;
  if (!ReactionsModifyTotal(p) && total.IsValid(majorant_temperature)) {
    return total.xs.at(std::get<ContinuousEnergy>(p.GetEnergy()));
  }
  else {
    return std::accumulate(
        reactions.cbegin(), reactions.cend(), MicroscopicCrossSection{0},
        [&p](const auto& accumulated, const auto& reaction) {
          return accumulated + reaction->GetMajorant(p);
        });
  }
}

MicroscopicCrossSection Continuous::GetTotal(const Particle& p) const noexcept {
  const auto requested_temperature =
      p.GetCell().temperature->at(p.GetPosition());
  if (!ReactionsModifyTotal(p) && total.IsValid(requested_temperature)) {
    return total.xs.at(std::get<ContinuousEnergy>(p.GetEnergy()));
  }
  else {
    return std::accumulate(
        reactions.cbegin(), reactions.cend(), MicroscopicCrossSection{0},
        [&p](const auto& accumulated, const auto& reaction) {
          return accumulated + reaction->GetCrossSection(p);
        });
  }
}

void Continuous::Interact(Particle& p) const noexcept {
  const MicroscopicCrossSection threshold =
      std::uniform_real_distribution{}(p.rng) * GetTotal(p);
  MicroscopicCrossSection accumulated{0};
  for (const auto& candidate : reactions) {
    accumulated += candidate->GetCrossSection(p);
    if (accumulated > threshold) {
      candidate->Interact(p);
      return;
    }
  }
  // If no reaction found, resample tail-recursively
  return Interact(p);
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

bool Continuous::ReactionsModifyTotal(const Particle& p) const noexcept {
  return std::any_of(
      reactions.cbegin(), reactions.cend(),
      [&p](const auto& reaction) { return reaction->ModifiesTotal(p); });
}
