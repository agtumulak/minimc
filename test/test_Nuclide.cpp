#include "NuclearData.hpp"
#include "Nuclide.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "Statistics.hpp"
#include "XMLDocument.hpp"
#include "World.hpp"
#include "catch2/catch.hpp"

#include <algorithm>
#include <random>
#include <variant>

TEST_CASE("nonexistent nuclide name throws exception") {
  XMLDocument doc{"simple_multigroup.xml"};
  REQUIRE_THROWS_WITH(
      Nuclide(doc.root, "nonexistent"),
      "Nuclide node \"nonexistent\" not found. Must be one of: [\"hydrogen\", "
      "\"oxygen\", \"uranium235\", \"inconsistent capture\", \"inconsistent "
      "scatter\", \"wrong particle\", ]");
}

TEST_CASE("poorly formed multigroup data for Nuclide throws exception") {
  XMLDocument doc{"simple_multigroup.xml"};
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
  SECTION("Multigroup methods") {
    XMLDocument doc{"simple_multigroup.xml"};
    World w{doc.root};
    Particle neutron_group1{
        Point{}, Direction{1, 0, 0}, Group{1}, Particle::Type::neutron};
    neutron_group1.SetCell(w.FindCellContaining(neutron_group1.GetPosition()));
    Particle neutron_group2{
        Point{}, Direction{1, 0, 0}, Group{2}, Particle::Type::neutron};
    neutron_group2.SetCell(w.FindCellContaining(neutron_group2.GetPosition()));
    const Nuclide hydrogen{doc.root, "hydrogen"};
    const Nuclide oxygen{doc.root, "oxygen"};
    const Nuclide uranium235{doc.root, "uranium235"};

    REQUIRE(hydrogen.GetTotal(neutron_group1) == 1);
    REQUIRE(hydrogen.GetTotal(neutron_group2) == 1);
    REQUIRE(oxygen.GetTotal(neutron_group1) == 1.5);
    REQUIRE(oxygen.GetTotal(neutron_group2) == 1);
    REQUIRE(uranium235.GetTotal(neutron_group1) == 1.33);
    REQUIRE(uranium235.GetTotal(neutron_group2) == 2.67);

    SECTION("SampleReaction() returns expected number of scatters") {
      // purely scattering for Group 1
      REQUIRE(
          ScatterProbability(rng, hydrogen, neutron_group1) ==
          Approx(1).epsilon(epsilon::Bernoulli(1, samples)));
      // purely absorbing for Group 2
      REQUIRE(
          ScatterProbability(rng, hydrogen, neutron_group2) ==
          Approx(0).epsilon(epsilon::Bernoulli(0, samples)));
      // total microscopic cross section is 1.5
      // scattering cross section is 1.0
      REQUIRE(
          ScatterProbability(rng, oxygen, neutron_group1) ==
          Approx(1. / 1.5).epsilon(epsilon::Bernoulli(1. / 1.5, samples)));
      // total microscopic cross section is 1.0
      // scattering cross section is 0.5
      REQUIRE(
          ScatterProbability(rng, oxygen, neutron_group2) ==
          Approx(0.5).epsilon(epsilon::Bernoulli(0.5, samples)));
      // total microscopic cross section is 1.33
      // scattering cross section is 1
      REQUIRE(
          ScatterProbability(rng, uranium235, neutron_group1) ==
          Approx(1. / 1.33).epsilon(epsilon::Bernoulli(1. / 1.33, samples)));
      // total microscopic cross section is 2.67
      // scattering cross section is 1
      REQUIRE(
          ScatterProbability(rng, uranium235, neutron_group2) ==
          Approx(1. / 2.67).epsilon(epsilon::Bernoulli(1. / 2.67, samples)));
    }
    SECTION("Scatter() returns expected number of Particles in lower energy "
            "Group") {
      auto LowerEnergyProbability = [](std::minstd_rand& rng, const Nuclide& n,
                                       const Particle& p) {
        size_t lower_energy{0};
        for (size_t i = 1; i <= samples; i++) {
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

      // The case LowerEnergyProbability(rng, hydrogen, neutron_group2) is not
      // tested. The particle neutron_group2 would never call Scatter() as it
      // has zero scatter cross section.

      // equally likely to inscatter (Group 1 -> Group 1) or downscatter
      // (Group 1 -> Group 2)
      REQUIRE(
          LowerEnergyProbability(rng, oxygen, neutron_group1) ==
          Approx(0.5).epsilon(epsilon::Bernoulli(0.5, samples)));
      // purely inscattering (Group 2 -> Group 2)
      REQUIRE(
          LowerEnergyProbability(rng, oxygen, neutron_group2) ==
          Approx(1).epsilon(epsilon::Bernoulli(1, samples)));
      // purely inscatteirng (Group 1 -> Group 1)
      REQUIRE(
          LowerEnergyProbability(rng, uranium235, neutron_group1) ==
          Approx(0).epsilon(epsilon::Bernoulli(0, samples)));
      // purely inscattering (Group 2 -> Group 2)
      REQUIRE(
          LowerEnergyProbability(rng, uranium235, neutron_group2) ==
          Approx(1).epsilon(epsilon::Bernoulli(1, samples)));
    }
    SECTION("Fission() returns expected number and spectrum of Particles") {
      struct FissionStatisticsResult {
        double ceil_nu_prob;
        double lower_energy_prob;
      };
      auto FissionStatistics = [](RNG& rng, const Nuclide& n,
                                  const Particle& p) {
        std::vector<Particle> bank;
        for (size_t i = 1; i <= samples; i++) {
          auto p_copy = p;
          const auto yield = n.Fission(rng, p_copy);
          for (const auto& fission_particle : n.Fission(rng, p_copy)) {
            bank.push_back(fission_particle);
          }
        }
        const auto sample_nu = static_cast<double>(bank.size()) / samples;
        const auto ceil_nu_prob =
            sample_nu - static_cast<unsigned int>(sample_nu);
        const auto lower_energy_prob =
            static_cast<double>(std::count_if(
                bank.cbegin(), bank.cend(),
                [](const auto& p) {
                  return std::get<Group>(p.GetEnergy()) == Group{2};
                })) /
            bank.size();
        return FissionStatisticsResult{ceil_nu_prob, lower_energy_prob};
      };

      // The case FissionStatistics(rng, uranium235, neutron_group1) is not
      // tested. The particle neutron_group1 would never call Fission() as it
      // has zero fission cros section.

      // Only neutron_group2 can cause a fission
      const auto stats = FissionStatistics(rng, uranium235, neutron_group2);
      // nu bar is 2.43
      REQUIRE(
          stats.ceil_nu_prob ==
          Approx(0.43).epsilon(epsilon::Bernoulli(0.43, samples)));
      // chi spectrum gives equal likelihood of fission particle in Group 1
      // or Group 2
      REQUIRE(
          stats.lower_energy_prob ==
          Approx(0.5).epsilon(epsilon::Bernoulli(0.5, samples)));
    }
  }
  SECTION("Continuous methods") {
    XMLDocument doc{"simple_continuous.xml"};
    World w{doc.root};
    Particle neutron{
        Point{}, Direction{1, 0, 0}, ContinuousEnergy{0.999e-6},
        Particle::Type::neutron};
    neutron.SetCell(w.FindCellContaining(neutron.GetPosition()));
    const Nuclide hydrogen{doc.root, "hydrogen"};
    const Nuclide oxygen{doc.root, "oxygen"};
    const Nuclide uranium235{doc.root, "uranium235"};
    REQUIRE_THROWS_WITH(
        Nuclide(doc.root, "badpath"), "File not found: /not/real");
    REQUIRE(hydrogen.GetTotal(neutron) == 20.74762);
    REQUIRE(oxygen.GetTotal(neutron) == 3.796956);
    REQUIRE(uranium235.GetTotal(neutron) == 92.60272);
    SECTION("SampleReaction() returns expected number of scatters") {
      // At 1.0 eV, for hydrogen,
      // capture: 0.05293893
      // elastic: 20.69468
      // total: 20.74762
      constexpr auto h_expected = 20.69468 / 20.74762;
      REQUIRE(
          ScatterProbability(rng, hydrogen, neutron) ==
          Approx(h_expected).epsilon(epsilon::Bernoulli(h_expected, samples)));
      // At 1.0 eV, for oxygen,
      // capture: 2.715658E-5
      // elastic: 3.796929
      // total: 3.796956
      constexpr auto o_expected = 3.796929 / 3.796956;
      REQUIRE(
          ScatterProbability(rng, oxygen, neutron) ==
          Approx(o_expected).epsilon(epsilon::Bernoulli(o_expected, samples)));
      // At 1.0 eV, for uranium-235,
      // capture: 12.51272
      // elastic: 12.63409
      // total: 92.60272
      constexpr auto u_expected = 12.63409 / 92.60272;
      REQUIRE(
          ScatterProbability(rng, uranium235, neutron) ==
          Approx(u_expected).epsilon(epsilon::Bernoulli(u_expected, samples)));
    }
    SECTION("Fission() returns expected number and spectrum of Particles") {
      struct FissionStatisticsResult {
        double ceil_nu_prob;
        double lower_energy_prob;
      };
      auto FissionStatistics = [](RNG& rng, const Nuclide& n, const Particle& p,
                                  const ContinuousEnergy& chi_median) {
        std::vector<Particle> bank;
        for (size_t i = 1; i <= samples; i++) {
          auto p_copy = p;
          for (const auto& fission_particle : n.Fission(rng, p_copy)) {
            bank.push_back(fission_particle);
          }
        }
        const auto sample_nu = static_cast<double>(bank.size()) / samples;
        const auto ceil_nu_prob =
            sample_nu - static_cast<unsigned int>(sample_nu);
        const auto lower_energy_prob =
            static_cast<double>(std::count_if(
                bank.cbegin(), bank.cend(),
                [&chi_median](const auto& p) {
                  return std::get<ContinuousEnergy>(p.GetEnergy()) <=
                         chi_median;
                })) /
            bank.size();
        return FissionStatisticsResult{ceil_nu_prob, lower_energy_prob};
      };
      // At 1.0 eV, for uranium235,
      // nu bar: 2.43385
      // median of chi: between 1.6 MeV and 1.7 MeV
      const auto chi_median = 1.65;
      const auto stats =
          FissionStatistics(rng, uranium235, neutron, chi_median);
      REQUIRE(
          stats.ceil_nu_prob ==
          Approx(0.43385).epsilon(epsilon::Bernoulli(0.43385, samples)));
      REQUIRE(
          stats.lower_energy_prob ==
          Approx(0.5).epsilon(epsilon::Bernoulli(0.5, samples)));
    }
  }
}
