#include "Source.hpp"

#include "pugixml.hpp"

#include <cassert>
#include <iosfwd>
#include <string>

// Template class specialization instantiations

template class ConstantDistribution<Point>;
template class ConstantDistribution<Direction>;
template class ConstantDistribution<Energy>;
template class ConstantDistribution<Particle>;

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
  else if constexpr (std::is_same_v<T, Particle>) {
    if (distribution_name == "constant") {
      const auto particle_type =
          ToParticle(distribution_node.attribute("type").as_string());
      return std::make_unique<ConstantDistribution<T>>(particle_type);
    }
  }
  else {
    constexpr bool is_distribution =
        std::is_same_v<T, Point> || std::is_same_v<T, Direction> ||
        std::is_same_v<T, Energy> || std::is_same_v<T, Particle>;
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

// Source

//// public

Source::Source(const pugi::xml_node& source_node)
    : position{Distribution<Point>::Create(source_node.child("position"))},
      direction{
          Distribution<Direction>::Create(source_node.child("direction"))},
      energy{Distribution<Energy>::Create(source_node.child("energy"))},
      particle{Distribution<Particle>::Create(source_node.child("particle"))} {}

Point Source::SamplePosition(RNG& rng) const noexcept {
  return position->Sample(rng);
}

Direction Source::SampleDirection(RNG& rng) const noexcept {
  return direction->Sample(rng);
}

Energy Source::SampleEnergy(RNG& rng) const noexcept {
  return energy->Sample(rng);
}

Particle Source::SampleParticle(RNG& rng) const noexcept {
  return particle->Sample(rng);
}
