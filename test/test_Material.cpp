#include "Material.hpp"
#include "Statistics.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

#include <random>

TEST_CASE("nonexistent material name throws exception") {
  XMLDocument doc{"simple.xml"};
  REQUIRE_THROWS_WITH(
      Material::FindNode(doc.root, "nonexistent"),
      Catch::Matchers::Equals(
          "Material node \"nonexistent\" not found. "
          "Must be one of: [\"water\", \"hydrogen\", \"oxygen\", ]"));
}

TEST_CASE("Material member methods work properly") {
  XMLDocument doc{"simple.xml"};
  const World w{doc.root};
  const auto& pit = w.cells.at(0);
  const auto& inner_shell = w.cells.at(1);
  const auto& hydrogen = *pit.material;
  const auto& water = *inner_shell.material;
  const Particle neutron_group1{Group{1}, Particle::Type::neutron};
  const Particle neutron_group2{Group{2}, Particle::Type::neutron};
  SECTION("SampleNuclide() returns expected values") {
    std::minstd_rand rng{42};
    const size_t samples{1000};
    auto HydrogenProbability = [&hydrogen](
                                   std::minstd_rand& rng, const Material& m,
                                   const Particle& p) {
      size_t hydrogen_count{0};
      for (size_t i = 1; i <= samples; i++) {
        if (m.SampleNuclide(rng, p).name == hydrogen.name) {
          hydrogen_count++;
        }
      }
      return static_cast<double>(hydrogen_count) / samples;
    };
    // hydrogen is the only nuclide in the hydrogen material
    REQUIRE(
        HydrogenProbability(rng, hydrogen, neutron_group1) ==
        Approx(1).epsilon(epsilon::Bernoulli(1, samples)));
    // hydrogen is the only nuclide in the hydrogen material
    REQUIRE(
        HydrogenProbability(rng, hydrogen, neutron_group2) ==
        Approx(1).epsilon(epsilon::Bernoulli(1, samples)));
    // total microscopic cross section is 1.165 for water; 0.67 for hydrogen
    // in water
    REQUIRE(
        HydrogenProbability(rng, water, neutron_group1) ==
        Approx(0.67 / 1.165)
            .epsilon(epsilon::Bernoulli(0.67 / 1.165, samples)));
    // total microscopic cross section is 1.0 for water; 0.67 for hydrogen in
    // water
    REQUIRE(
        HydrogenProbability(rng, water, neutron_group2) ==
        Approx(0.67).epsilon(epsilon::Bernoulli(0.67, samples)));
  }
}
