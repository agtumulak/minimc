#include "InteractionDelegate.hpp"

#include "ContinuousReaction.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "Reaction.hpp"
#include "pugixml.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <variant>

// InteractionDelegate

//// public

InteractionDelegate::~InteractionDelegate() noexcept {}

// Continuous

//// public

Continuous::Continuous(
    const pugi::xml_node& particle_node, const Nuclide& target)
    : // IIFE
      tnsl{[&particle_node, &target]() -> std::optional<ThermalScattering> {
        const auto& tnsl_node = particle_node.child("scatter").child("tnsl");
        if (!tnsl_node) {
          return std::nullopt;
        }
        if (Particle::ToType(particle_node.name()) != Particle::Type::neutron) {
          throw std::runtime_error(
              tnsl_node.path() +
              ": Only neutrons may have a thermal scattering library node");
        }
        return ThermalScattering{tnsl_node, target};
      }()},
      // IIFE
      reactions{[this, &particle_node, &target]() {
        std::vector<std::unique_ptr<const ContinuousReaction>> reactions;
        for (const auto& reaction_node : particle_node) {
          const std::string reaction_name = reaction_node.name();
          reactions.emplace_back(
              ContinuousReaction::Create(reaction_node, target, tnsl));
        }
        return reactions;
      }()} {}

const std::optional<ThermalScattering>& Continuous::GetTNSL() const {
  return tnsl;
}

MicroscopicCrossSection
Continuous::GetCellMajorant(const Particle& p) const noexcept {
  return std::accumulate(
      reactions.cbegin(), reactions.cend(), MicroscopicCrossSection{0},
      [&p](const auto& accumulated, const auto& reaction) {
        return accumulated + reaction->GetCellMajorant(p);
      });
}

void Continuous::Interact(
    Particle& p,
    std::vector<Estimator::Proxy>& estimator_proxies) const noexcept {
  const MicroscopicCrossSection threshold =
      std::uniform_real_distribution{}(p.rng) * GetCellMajorant(p);
  MicroscopicCrossSection accumulated{0};
  for (const auto& candidate : reactions) {
    accumulated += candidate->GetCellMajorant(p);
    if (accumulated > threshold) {
      candidate->Interact(p, estimator_proxies);
      return;
    }
  }
  // If no reaction found, resample tail-recursively
  return Interact(p, estimator_proxies);
}

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

const std::optional<ThermalScattering>& Multigroup::GetTNSL() const {
  throw std::runtime_error(
      "Thermal neutron scattering law does not exist for multigroup data");
}

MicroscopicCrossSection
Multigroup::GetCellMajorant(const Particle& p) const noexcept {
  return total.at(std::get<Group>(p.GetEnergy()));
}

void Multigroup::Interact(
    Particle& p,
    std::vector<Estimator::Proxy>& estimator_proxies) const noexcept {
  const MicroscopicCrossSection threshold{
      std::uniform_real_distribution{}(p.rng) * GetCellMajorant(p)};
  MicroscopicCrossSection accumulated{0};
  for (const auto& [reaction, xs] : reactions) {
    accumulated += xs.at(std::get<Group>(p.GetEnergy()));
    if (accumulated > threshold) {
      switch (reaction) {
      case Reaction::birth:
        assert(false); // Interact() only called after Stream()
        break;
      case Reaction::capture:
        Capture(p, estimator_proxies);
        break;
      case Reaction::scatter:
        Scatter(p, estimator_proxies);
        break;
      case Reaction::fission:
        Fission(p, estimator_proxies);
        break;
      case Reaction::leak:
        assert(false); // Interact() should never be called on leaked Particle
        break;
      }
      return;
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
  std::vector<MicroscopicCrossSection> column_sums;
  auto scatter_matrix = TwoDimensional{scatter_node};
  for (const auto& column : scatter_matrix) {
    column_sums.push_back(std::accumulate(
        column.begin(), column.end(), MicroscopicCrossSection{0}));
  }
  return column_sums;
}

Multigroup::ReactionsMap
Multigroup::CreateReactions(const pugi::xml_node& particle_node) {
  Multigroup::ReactionsMap reactions;
  for (const auto& reaction_node : particle_node) {
    switch (const auto reaction = ToReaction(reaction_node.name())) {
    case Reaction::birth:
      assert(false); // Interact() only called after StreamToNextCollision()
      break;
    case Reaction::capture:
      reactions.emplace(reaction, reaction_node);
      break;
    case Reaction::scatter:
      reactions.emplace(reaction, CreateScatterXS(reaction_node));
      break;
    case Reaction::fission:
      // fission cross sections "xs" are bundled along with "chi" and "nubar"
      reactions.emplace(reaction, reaction_node.child("xs"));
      break;
    case Reaction::leak:
      assert(false); // Interact() should never be called on leaked Particle
      break;
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

void Multigroup::Capture(
    Particle& p, std::vector<Estimator::Proxy>&) const noexcept {
  p.reaction = Reaction::capture;
}

void Multigroup::Scatter(
    Particle& p, std::vector<Estimator::Proxy>&) const noexcept {
  p.reaction = Reaction::scatter;
  const Real threshold = std::uniform_real_distribution{}(p.rng);
  Real accumulated{0};
  const auto& group_probs{
      scatter_probs.value().at(std::get<Group>(p.GetEnergy()))};
  for (Group g = 1; g <= max_group; g++) {
    accumulated += group_probs.at(g);
    if (accumulated > threshold) {
      p.SetEnergy(g);
      p.SetDirectionIsotropic();
      return;
    }
  }
  // scatter_probs is a NormalizedTwoDimensional so this shouldn't ever happen
  assert(false);
}

void Multigroup::Fission(
    Particle& p, std::vector<Estimator::Proxy>&) const noexcept {
  p.reaction = Reaction::fission;
  auto incident_group{std::get<Group>(p.GetEnergy())};
  // rely on the fact that double to int conversions essentially do a floor()
  size_t fission_yield(
      nubar.value().at(incident_group) +
      std::uniform_real_distribution{}(p.rng));
  const auto& group_probs{chi.value().at(incident_group)};
  for (size_t i = 0; i < fission_yield; i++) {
    const Real threshold = std::uniform_real_distribution{}(p.rng);
    Real accumulated{0};

    for (Group g = 1; g <= max_group; g++) {
      accumulated += group_probs.at(g);
      if (accumulated > threshold) {
        p.BankSecondaries(Direction{p.rng}, g);
        break;
      }
    }
  }
}
