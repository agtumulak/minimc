#include "Multigroup.hpp"

#include "Constants.hpp"

#include <algorithm>
#include <functional>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <variant>

// Multigroup

//// public

Multigroup::Multigroup(const pugi::xml_node& particle_node)
    : nubar{particle_node.child("fission").child("nubar")
              ? std::make_optional<OneDimensional>(
                  particle_node.child("fission").child("nubar"))
              : std::nullopt},
      chi{particle_node.child("fission").child("chi")
              ? std::make_optional<NormalizedTwoDimensional>(
                  particle_node.child("fission").child("chi"))
              : std::nullopt},
      scatter_probs{particle_node.child("scatter")
              ? std::make_optional<NormalizedTwoDimensional>(
                  particle_node.child("scatter"))
              : std::nullopt},
      reactions{CreateReactions(particle_node)},
      total{CreateTotalXS(particle_node, reactions)},
      max_group{GroupStructureSize(particle_node.root())} {}

NuclearData::CrossSection
Multigroup::GetTotal(const Particle& p) const noexcept {
  return total.at(std::get<Group>(p.GetEnergy()));
}

NuclearData::CrossSection
Multigroup::GetFission(const Particle& p) const noexcept {
  try {
    return reactions.at(NuclearData::Reaction::fission)
        .at(std::get<Group>(p.GetEnergy()));
  }
  catch (const std::out_of_range& e) {
    return 0;
  }
}

Real Multigroup::GetNuBar(const Particle& p) const noexcept {
  if (nubar) {
    return nubar->at(std::get<Group>(p.GetEnergy()));
  }
  else {
    return 0;
  }
}

void Multigroup::Scatter(RNG& rng, Particle& p) const noexcept {
  const Real threshold = std::uniform_real_distribution{}(rng);
  Real accumulated{0};
  const auto& group_probs{
      scatter_probs.value().at(std::get<Group>(p.GetEnergy()))};
  for (Group g = 1; g <= max_group; g++) {
    accumulated += group_probs.at(g);
    if (accumulated > threshold) {
      p.SetEnergy(g);
      p.SetDirectionIsotropic(rng);
      return;
    }
  }
  // scatter_probs is a NormalizedTwoDimensional so this shouldn't ever happen
  assert(false);
}

std::vector<Particle>
Multigroup::Fission(RNG& rng, Particle& p) const noexcept {
  std::vector<Particle> fission_neutrons;
  p.Kill();
  auto incident_group{std::get<Group>(p.GetEnergy())};
  // rely on the fact that double to int conversions essentially do a floor()
  size_t fission_yield(
      nubar.value().at(incident_group) + std::uniform_real_distribution{}(rng));
  const auto& group_probs{chi.value().at(incident_group)};
  for (size_t i = 0; i < fission_yield; i++) {
    const Real threshold = std::uniform_real_distribution{}(rng);
    Real accumulated{0};
    for (Group g = 1; g <= max_group; g++) {
      accumulated += group_probs.at(g);
      if (accumulated > threshold) {
        fission_neutrons.emplace_back(
            p.GetPosition(), Direction::CreateIsotropic(rng), Energy{Group{g}},
            Particle::Type::neutron);
        fission_neutrons.back().SetCell(p.GetCell());
        // Each secondary produced by this fission must have a completely
        // deterministic history. Setting their seed here accomplishes this.
        fission_neutrons.back().seed =
            std::uniform_int_distribution<RNG::result_type>{1}(rng);
        break;
      }
    }
  }
  assert(fission_neutrons.size() == fission_yield);
  return fission_neutrons;
}

NuclearData::Reaction
Multigroup::SampleReaction(RNG& rng, const Particle& p) const noexcept {
  const CrossSection threshold{
      std::uniform_real_distribution{}(rng)*GetTotal(p)};
  CrossSection accumulated{0};
  for (const auto& [reaction, reaction_xs] : reactions) {
    accumulated += reaction_xs.at(std::get<Group>(p.GetEnergy()));
    if (accumulated > threshold) {
      return reaction;
    }
  }
  // TODO: Check each Nuclide has nonzero total cross section so this never
  //       happens
  assert(false);
}

// Multigroup::OneDimensional

//// public

Multigroup::OneDimensional::OneDimensional(const elements_type& elements)
    : elements{elements} {}

Multigroup::OneDimensional::OneDimensional(const pugi::xml_node& groupxs_node) {
  const Group G{GroupStructureSize(groupxs_node.root())};
  std::stringstream element_stringlist{groupxs_node.child_value()};
  Real element;
  while (element_stringlist >> element) {
    elements.push_back(element);
  }
  if (elements.size() != G) {
    throw std::runtime_error(
        groupxs_node.path() + ": Expected " + std::to_string(G) +
        " entries but got " + std::to_string(elements.size()));
  }
}

const Real& Multigroup::OneDimensional::at(const Group g) const {
  return elements.at(g - 1);
}

Multigroup::OneDimensional::elements_type::iterator
Multigroup::OneDimensional::begin() noexcept {
  return elements.begin();
}

Multigroup::OneDimensional::elements_type::const_iterator
Multigroup::OneDimensional::begin() const noexcept {
  return elements.begin();
}

Multigroup::OneDimensional::elements_type::iterator
Multigroup::OneDimensional::end() noexcept {
  return elements.end();
};

Multigroup::OneDimensional::elements_type::const_iterator
Multigroup::OneDimensional::end() const noexcept {
  return elements.end();
};

// Multigroup::TwoDimensional

//// public

Multigroup::TwoDimensional::TwoDimensional(const pugi::xml_node& groupxs_node) {
  const Group G{GroupStructureSize(groupxs_node.root())};
  // Load scatter matrix entries as flattened values
  std::vector<Real> flattened;
  std::stringstream element_stringlist{groupxs_node.child_value()};
  Real element;
  while (element_stringlist >> element) {
    flattened.push_back(element);
  }
  const size_t expected_size{G * G};
  if (flattened.size() != expected_size) {
    throw std::runtime_error(
        groupxs_node.path() + ": Expected " + std::to_string(expected_size) +
        " entries but got " + std::to_string(flattened.size()));
  }
  // Construct column by column
  for (Group gp = 1; gp <= G; gp++) {
    std::vector<Real> column;
    for (Group g = 1; g <= G; g++) {
      column.push_back(flattened.at(G * (g - 1) + (gp - 1)));
    }
    elements.push_back(column);
  }
}

const Multigroup::OneDimensional&
Multigroup::TwoDimensional::at(Group g) const {
  return elements.at(g - 1);
}

Multigroup::TwoDimensional::elements_type::const_iterator
Multigroup::TwoDimensional::begin() const noexcept {
  return elements.begin();
}

Multigroup::TwoDimensional::elements_type::const_iterator
Multigroup::TwoDimensional::end() const noexcept {
  return elements.end();
}

// Multigroup::NormalizedTwoDimensional

//// public

Multigroup::NormalizedTwoDimensional::NormalizedTwoDimensional(
    const pugi::xml_node& groupxs_node)
    : TwoDimensional(groupxs_node) {
  std::for_each(
      elements.begin(), elements.end(),
      [](elements_type::value_type& one_dimensional) {
        auto column_sum = std::accumulate(
            one_dimensional.begin(), one_dimensional.end(), Real{0},
            std::plus<Real>{});
        if (column_sum == 0.) {
          return;
        }
        else {
          std::for_each(
              one_dimensional.begin(), one_dimensional.end(),
              [&column_sum](auto& element) { element /= column_sum; });
        }
      });
}

// Multigroup

//// private

Group Multigroup::GroupStructureSize(const pugi::xml_node& root) noexcept {
  return root.root()
      .child("minimc")
      .child("nuclides")
      .child("multigroup")
      .attribute("groups")
      .as_uint();
}

Multigroup::OneDimensional
Multigroup::CreateScatterXS(const pugi::xml_node& scatter_node) {
  std::vector<CrossSection> column_sums;
  auto scatter_matrix = TwoDimensional{scatter_node};
  for (const auto& column : scatter_matrix) {
    column_sums.push_back(
        std::accumulate(column.begin(), column.end(), CrossSection{0}));
  }
  return column_sums;
}

Multigroup::ReactionsMap
Multigroup::CreateReactions(const pugi::xml_node& particle_node) {
  Multigroup::ReactionsMap reactions;
  for (const auto& reaction_node : particle_node) {
    const auto reaction{NuclearData::ToReaction(reaction_node.name())};
    if (reaction == NuclearData::Reaction::scatter) {
      reactions.emplace(reaction, CreateScatterXS(reaction_node));
    }
    else if (reaction == NuclearData::Reaction::fission){
      // fission cross sections "xs" are bundled along with "chi" and "nubar"
      reactions.emplace(reaction, reaction_node.child("xs"));
    }
    else {
      reactions.emplace(reaction, reaction_node);
    }
  }
  return reactions;
}

Multigroup::OneDimensional Multigroup::CreateTotalXS(
    const pugi::xml_node& particle_node,
    const ReactionsMap& reactions) noexcept {
  const Group G{GroupStructureSize(particle_node.root())};
  return std::accumulate(
      reactions.begin(), reactions.end(), std::vector<Real>(G, 0),
      [](auto& accumulated, const auto& reaction_xs) noexcept {
        std::transform(
            accumulated.begin(), accumulated.end(), reaction_xs.second.begin(),
            accumulated.begin(), std::plus<Real>{});
        return accumulated;
      });
}
