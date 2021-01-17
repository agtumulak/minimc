#include "Multigroup.hpp"

#include <algorithm>
#include <functional>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <variant>

// Multigroup

//// public

Multigroup::Multigroup(const pugi::xml_node& particle_node, const Group G)
    : scatter_matrix{CreateScatterMatrix(particle_node, G)},
      reactions{CreateReactions(particle_node, G, scatter_matrix)},
      total{CreateTotalXS(reactions, G)} {}

NuclearData::CrossSection
Multigroup::GetTotal(const Particle& p) const noexcept {
  return total.at(std::get<Group>(p.GetEnergy()));
}

void Multigroup::Scatter(std::minstd_rand& rng, Particle& p) const {
  const CrossSection threshold = std::uniform_real_distribution{}(
      rng)*GetReaction(p, NuclearData::Reaction::scatter);
  CrossSection accumulated{0};
  Group g{0};
  for (const auto& outgoing_probability : GetOutgoingScatterProbs(p)) {
    accumulated += outgoing_probability;
    g++;
    if (accumulated > threshold) {
      p.SetEnergy(g);
      p.SetDirectionIsotropic(rng);
      return;
    }
  }
  throw std::runtime_error(
      "SampleOutgoing reached end of possible outgoing Groups");
}

NuclearData::Reaction Multigroup::SampleReaction(
    std::minstd_rand& rng, const Particle& p) const noexcept {
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

Real& Multigroup::OneDimensional::at(const Group g) {
  return elements.at(g - 1);
}

const Real& Multigroup::OneDimensional::at(const Group g) const {
  return elements.at(g - 1);
}

Multigroup::OneDimensional::elements_type::const_iterator
Multigroup::OneDimensional::begin() const noexcept {
  return elements.begin();
}
Multigroup::OneDimensional::elements_type::const_iterator
Multigroup::OneDimensional::end() const noexcept {
  return elements.end();
};

// Multigroup::TwoDimensional

//// public

Multigroup::TwoDimensional::TwoDimensional(const elements_type& elements)
    : elements{elements} {}

Multigroup::OneDimensional& Multigroup::TwoDimensional::at(Group g) {
  return elements.at(g - 1);
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

// Multigroup

//// private

Multigroup::TwoDimensional Multigroup::CreateScatterMatrix(
    const pugi::xml_node& particle_node, const Group G) {
  TwoDimensional scatter_matrix(
      std::vector<OneDimensional>(G, std::vector<Real>(G, 0)));
  const auto& scatter_node{particle_node.child("scatter")};
  // Load scatter matrix entries as flattened values
  std::vector<Real> flattened;
  std::stringstream elements{scatter_node.child_value()};
  Real element;
  while (elements >> element) {
    flattened.push_back(element);
  }
  const size_t expected_size{G * G};
  if (flattened.size() != expected_size) {
    throw std::runtime_error(
        scatter_node.path() + ": Expected " + std::to_string(expected_size) +
        " entries but got " + std::to_string(flattened.size()));
  }
  // Unflatten matrix entries
  for (Group g = 1; g <= G; g++) {      // g corresponds to outgoing group
    for (Group gp = 1; gp <= G; gp++) { // gp corresponds to incoming group
      scatter_matrix.at(gp).at(g) = flattened.at(G * (g - 1) + (gp - 1));
    }
  }
  return scatter_matrix;
}

Multigroup::ReactionsMap Multigroup::CreateReactions(
    const pugi::xml_node& particle_node, const Group G,
    const TwoDimensional& scatter_matrix) {
  Multigroup::ReactionsMap reactions;
  for (const auto& reaction_node : particle_node) {
    std::vector<Real> xs;
    const auto reaction{NuclearData::ToReaction(reaction_node.name())};
    // scatter cross section is obtained by summing columns of scatter_matrix
    if (reaction == NuclearData::Reaction::scatter) {
      for (const auto& column : scatter_matrix) {
        xs.push_back(std::accumulate(
            column.begin(), column.end(), Real{0}, std::plus<Real>()));
      }
    }
    else {
      // everything else is a 1d vector
      std::stringstream elements{reaction_node.child_value()};
      Real element;
      while (elements >> element) {
        xs.push_back(element);
      }
      if (xs.size() != G) {
        throw std::runtime_error(
            reaction_node.path() + ": Expected " + std::to_string(G) +
            " entries but got " + std::to_string(xs.size()));
      }
    }
    reactions.emplace(reaction, xs);
  }
  return reactions;
}

Multigroup::OneDimensional Multigroup::CreateTotalXS(
    const ReactionsMap& reactions, const Group G) noexcept {
  return std::accumulate(
      reactions.begin(), reactions.end(), std::vector<Real>(G, 0),
      [](auto& accumulated, const auto& reaction_xs) noexcept {
        std::transform(
            accumulated.begin(), accumulated.end(), reaction_xs.second.begin(),
            accumulated.begin(), std::plus<Real>());
        return accumulated;
      });
}

NuclearData::CrossSection
Multigroup::GetReaction(const Particle& p, const Reaction r) const noexcept {
  return reactions.at(r).at(std::get<Group>(p.GetEnergy()));
}

const Multigroup::OneDimensional&
Multigroup::GetOutgoingScatterProbs(const Particle& p) const noexcept {
  return scatter_matrix.at(std::get<Group>(p.GetEnergy()));
}
