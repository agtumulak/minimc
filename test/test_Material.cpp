#include "Material.hpp"
#include "Statistics.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

#include <random>

TEST_CASE("nonexistent material name throws exception") {
  XMLDocument doc{"simple_multigroup.xml"};
  REQUIRE_THROWS_WITH(
      Material::FindNode(doc.root, "nonexistent"),
      Catch::Matchers::Equals(
          "Material node \"nonexistent\" not found. "
          "Must be one of: [\"water\", \"hydrogen\", \"oxygen\", ]"));
}

TEST_CASE("Material member methods work properly") {
  std::minstd_rand rng{42};
  const size_t samples{1000};
  auto HydrogenProbability = [](std::minstd_rand& rng, const Material& m,
                                const Particle& p) {
    size_t hydrogen_count{0};
    for (size_t i = 1; i <= samples; i++) {
      if (m.SampleNuclide(rng, p).name == "hydrogen") {
        hydrogen_count++;
      }
    }
    return static_cast<double>(hydrogen_count) / samples;
  };
  SECTION("SampleNuclide() returns expected values (multigroup)") {
    XMLDocument doc{"simple_multigroup.xml"};
    const World w{doc.root};
    const auto& hydrogen = *w.cells.at(0).material;
    const auto& water = *w.cells.at(1).material;
    const Particle neutron_group1{
        Point{}, Direction{1, 0, 0}, Group{1}, Particle::Type::neutron};
    const Particle neutron_group2{
        Point{}, Direction{1, 0, 0}, Group{2}, Particle::Type::neutron};
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
  SECTION("SampleNuclide() returns expected values (continuous)") {
    XMLDocument doc{"simple_continuous.xml"};
    const World w{doc.root};
    const auto& hydrogen = *w.cells.at(0).material;
    const auto& water = *w.cells.at(1).material;
    const Particle neutron{
        Point{}, Direction{1, 0, 0}, ContinuousEnergy{0.999},
        Particle::Type::neutron};
    // hydrogen is the only nuclide in the hydrogen material
    REQUIRE(
        HydrogenProbability(rng, hydrogen, neutron) ==
        Approx(1).epsilon(epsilon::Bernoulli(1, samples)));
    // At 1.0 MeV,
    // oxygen in water: afrac * sigma_t = 0.33 * 3.796956 = 1.25299548
    // hydrogen in water: afrac * sigma_t = 0.67 * 20.74762 = 13.9009054
    constexpr auto expected = 13.9009054 / (1.25299548 + 13.9009054);
    REQUIRE(
        HydrogenProbability(rng, water, neutron) ==
        Approx(expected).epsilon(epsilon::Bernoulli(expected, samples)));
  }
}
