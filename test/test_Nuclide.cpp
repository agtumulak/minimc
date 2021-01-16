#include "NuclearData.hpp"
#include "Nuclide.hpp"
#include "Particle.hpp"
#include "Statistics.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

#include <random>
#include <variant>

TEST_CASE("nonexistent nuclide name throws exception") {
  XMLDocument doc{"simple.xml"};
  REQUIRE_THROWS_WITH(
      Nuclide(doc.root, "nonexistent"),
      "Nuclide node \"nonexistent\" not found. Must be one of: [\"hydrogen\", "
      "\"oxygen\", \"inconsistent capture\", \"inconsistent scatter\", \"wrong "
      "particle\", ]");
}

TEST_CASE("poorly formed multigroup data for Nuclide throws exception") {
  XMLDocument doc{"simple.xml"};
  REQUIRE_THROWS_WITH(
      Nuclide(doc.root, "inconsistent capture"),
      "/minimc/nuclides/multigroup/nuclide/neutron/capture: Expected 2 entries "
      "but got 3");
  REQUIRE_THROWS_WITH(
      Nuclide(doc.root, "inconsistent scatter"),
      "/minimc/nuclides/multigroup/nuclide/neutron/scatter: Expected 4 entries "
      "but got 5");
  REQUIRE_THROWS_WITH(
      Nuclide(doc.root, "wrong particle"),
      "/minimc/nuclides/multigroup/nuclide: \"neutron\" node not found");
}

TEST_CASE("Nuclide member methods work properly") {
  XMLDocument doc{"simple.xml"};
  const Particle neutron_group1{Group{1}, Particle::Type::neutron};
  const Particle neutron_group2{Group{2}, Particle::Type::neutron};
  const Nuclide hydrogen{doc.root, "hydrogen"};
  const Nuclide oxygen{doc.root, "oxygen"};

  REQUIRE(hydrogen.GetTotal(neutron_group1) == 1);
  REQUIRE(hydrogen.GetTotal(neutron_group2) == 1);
  REQUIRE(oxygen.GetTotal(neutron_group1) == 1.5);
  REQUIRE(oxygen.GetTotal(neutron_group2) == 1);

  SECTION("SampleReaction() returns expected number of scatters") {
    std::minstd_rand rng{42};
    const size_t samples{1000};
    auto ScatterProbability = [](std::minstd_rand& rng, const Nuclide& n,
                                 const Particle& p) {
      size_t scatters{0};
      for (size_t i = 1; i <= samples; i++) {
        if (n.SampleReaction(rng, p) == NuclearData::Reaction::scatter) {
          scatters++;
        }
      };
      return static_cast<double>(scatters) / samples;
    };
    // purely scattering for Group 1
    REQUIRE(
        ScatterProbability(rng, hydrogen, neutron_group1) ==
        Approx(1).epsilon(epsilon::Bernoulli(1, samples)));
    // purely absorbing for Group 2
    REQUIRE(
        ScatterProbability(rng, hydrogen, neutron_group2) ==
        Approx(0).epsilon(epsilon::Bernoulli(0, samples)));
    // total microscopic cross section is 1.5 for Group 1
    // scattering cross section is 1.0 for Group 1
    REQUIRE(
        ScatterProbability(rng, oxygen, neutron_group1) ==
        Approx(1. / 1.5).epsilon(epsilon::Bernoulli(1. / 1.5, samples)));
    // total microscopic cross section is 1.0 for Group 2
    // scattering cross section is 0.5 for Group 2
    REQUIRE(
        ScatterProbability(rng, oxygen, neutron_group2) ==
        Approx(0.5).epsilon(epsilon::Bernoulli(0.5, samples)));
  }
  SECTION(
      "Scatter() returns expected number of Particles in lower energy Group") {
    std::minstd_rand rng{42};
    const size_t samples{1000};
    auto LowerEnergyProbability = [](std::minstd_rand& rng, const Nuclide& n,
                                     const Particle& p) {
      size_t lower_energy{0};
      for (size_t i = 1; i <= 1000; i++) {
        auto p_copy = p;
        n.Scatter(rng, p_copy);
        if (std::get<Group>(p_copy.GetEnergy()) == Group{2}) {
          lower_energy++;
        }
      };
      return static_cast<double>(lower_energy) / samples;
    };
    // purely downscattering (Group 1 -> Group 2)
    REQUIRE(
        LowerEnergyProbability(rng, hydrogen, neutron_group1) ==
        Approx(1).epsilon(epsilon::Bernoulli(1, samples)));
    // zero scattering cross section; shouldn't be called during transport
    REQUIRE_THROWS_WITH(
        LowerEnergyProbability(rng, hydrogen, neutron_group2),
        "SampleOutgoing reached end of possible outgoing Groups");
    // equally likely to inscatter (Group 1 -> Group 1) or downscatter
    // (Group 1 -> Group 2)
    REQUIRE(
        LowerEnergyProbability(rng, oxygen, neutron_group1) ==
        Approx(0.5).epsilon(epsilon::Bernoulli(0.5, samples)));
    // purely inscattering (Group 2 -> Group 2)
    REQUIRE(
        LowerEnergyProbability(rng, oxygen, neutron_group2) ==
        Approx(1).epsilon(epsilon::Bernoulli(1, samples)));
  }
}
