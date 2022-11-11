#include "Source.hpp"

#include "BasicTypes.hpp"
#include "Constants.hpp"
#include "World.hpp"
#include "pugixml.hpp"

#include <cassert>
#include <cmath>
#include <iosfwd>
#include <random>
#include <string>
#include <type_traits>
#include <unordered_map>

// Template class specialization instantiations

template class ConstantDistribution<Point>;
template class ConstantDistribution<Direction>;
template class ConstantDistribution<Energy>;
template class ConstantDistribution<Particle::Type>;

// Distribution

//// public

template <typename T>
std::unique_ptr<const Distribution<T>>
Distribution<T>::Create(const pugi::xml_node& property_node) {
  // any assert(false) statements mean that the validator should have caught it
  const auto& distribution_node = property_node.first_child();
  const std::string distribution_name = distribution_node.name();
  // compile-time type checks wow this is truly the future
  if constexpr (std::is_same_v<T, Point>) {
    if (distribution_name == "constant") {
      return std::make_unique<ConstantDistribution<T>>(T{distribution_node});
    }
    else {
      assert(false);
    }
  }
  else if constexpr (std::is_same_v<T, Direction>) {
    if (distribution_name == "constant") {
      return std::make_unique<ConstantDistribution<T>>(T{distribution_node});
    }
    else if (distribution_name == "isotropic") {
      return std::make_unique<IsotropicDistribution>();
    }
    else if (distribution_name == "isotropic-flux") {
      return std::make_unique<IsotropicFlux>(distribution_node);
    }
    else {
      assert(false);
    }
  }
  else if constexpr (std::is_same_v<T, Energy>) {
    if (distribution_name == "constant") {
      const std::string energy_type = property_node.root()
                                          .child("minimc")
                                          .child("nuclides")
                                          .first_child()
                                          .name();
      const std::string energy_value =
          distribution_node.attribute("energy").as_string();
      if (energy_type == "multigroup") {
        return std::make_unique<ConstantDistribution<T>>(
            Group{std::stoull(energy_value)});
      }
      else if (energy_type == "continuous") {
        return std::make_unique<ConstantDistribution<T>>(
            ContinuousEnergy{std::stod(energy_value)});
      }
      else {
        assert(false);
      }
    }
  }
  else if constexpr (std::is_same_v<T, Particle::Type>) {
    if (distribution_name == "constant") {
      const auto particle_type =
          Particle::ToType(distribution_node.attribute("type").as_string());
      return std::make_unique<ConstantDistribution<T>>(particle_type);
    }
  }
  else {
    constexpr bool is_distribution =
        std::is_same_v<T, Point> || std::is_same_v<T, Direction> ||
        std::is_same_v<T, Energy> || std::is_same_v<T, Particle::Type>;
    static_assert(is_distribution, "Template parameter not supported");
  }
  assert(false);
  return {};
}

template <typename T> Distribution<T>::~Distribution() noexcept {};

// ConstantDistribution

//// public

template <typename T>
ConstantDistribution<T>::ConstantDistribution(const T& constant)
    : constant{constant} {};

template <typename T> T ConstantDistribution<T>::Sample(RNG&) const noexcept {
  return constant;
};

// IsotropicDistribution

//// public

Direction IsotropicDistribution::Sample(RNG& rng) const noexcept {
  return Direction{rng};
}

// IsotropicFlux

//// public

IsotropicFlux::IsotropicFlux(const pugi::xml_node& isotropic_flux_node) noexcept
    : reference{isotropic_flux_node} {}

Direction IsotropicFlux::Sample(RNG& rng) const noexcept {
  // sample a value of mu in [0, 1) with probability p(mu) = 2 * mu
  const Real mu = std::sqrt(std::uniform_real_distribution{}(rng));
  const Real phi = std::uniform_real_distribution{0., 2 * constants::pi}(rng);
  return Direction{reference, mu, phi};
}

// Source

//// public

Source::Source(
    const pugi::xml_node& source_node, const World& world,
    const std::vector<std::unique_ptr<const Perturbation>>& perturbations)
    : position{Distribution<Point>::Create(source_node.child("position"))},
      direction{
          Distribution<Direction>::Create(source_node.child("direction"))},
      energy{Distribution<Energy>::Create(source_node.child("energy"))},
      particle_type{Distribution<Particle::Type>::Create(
          source_node.child("particletype"))},
      world{world}, perturbations{perturbations} {}

Particle Source::Sample(RNG::result_type seed) const noexcept {
  // initialize indirect effects for new Particle
  std::unordered_map<const Perturbation*, Real> indirect_effects;
  for (const auto& perturbation : perturbations) {
    indirect_effects[perturbation.get()] = 0;
  }
  // evaluation order of arguments is undefined so do evaluation here
  RNG rng{seed};
  const auto sampled_position = position->Sample(rng);
  const auto sampled_direction = direction->Sample(rng);
  const auto sampled_energy = energy->Sample(rng);
  const auto sampled_particle_type = particle_type->Sample(rng);
  const auto sampled_seed = rng();
  return Particle{
      indirect_effects,
      sampled_position,
      sampled_direction,
      sampled_energy,
      &world.FindCellContaining(sampled_position),
      sampled_particle_type,
      sampled_seed};
}
