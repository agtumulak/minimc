#include "Cell.hpp"
#include "Statistics.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

#include <random>

TEST_CASE("Cell is constructed properly") {
  XMLDocument doc{"simple.xml"};
  const World w{doc.root};

  const auto& pit = w.cells.at(0);
  const auto& inner_shell = w.cells.at(1);
  const Particle neutron_group1{Group{1}, Particle::Type::neutron};
  const Particle neutron_group2{Group{2}, Particle::Type::neutron};

  SECTION("Cell returns correct nearest CSGSurface"){
    const auto& [surface, distance] =
        inner_shell.NearestSurface(Point{1.5, 0, 0}, Point{1, 0, 0});
    REQUIRE(surface->name == "middle sphere");
    REQUIRE(distance == Approx(0.5));
  }

  SECTION("Sampled distance to collision is close to expected value") {
    // define helper function for calling Cell::SampleDistance
    std::minstd_rand rng{42};
    const size_t samples{1000};
    auto TestSampleDistance =
        [&samples](std::minstd_rand& rng, const Cell& c, const Particle& p) {
          Real accumulated{0};
          for (size_t i = 1; i <= samples; i++) {
            accumulated += c.SampleCollisionDistance(rng, p);
          }
          return accumulated / samples;
        };
    // Total macroscopic cross section in pit is 1 for Group 1
    REQUIRE(
        TestSampleDistance(rng, pit, neutron_group1) ==
        Approx(1).epsilon(epsilon::Exponential(samples)));
    // Total macroscopic cross section in pit is 1 for Group 2
    REQUIRE(
        TestSampleDistance(rng, pit, neutron_group2) ==
        Approx(1).epsilon(epsilon::Exponential(samples)));
    // Total cross section in inner_shell for Group 1 is
    // 2 * (0.67 * 1 + 0.33 * 1.5) = 2.33
    REQUIRE(
        TestSampleDistance(rng, inner_shell, neutron_group1) ==
        Approx(1. / 2.33).epsilon(epsilon::Exponential(samples)));
    // Total cross section in inner_shell for Group 2 is
    // 2 * (0.67 * 1 + 0.33 * 1) = 2
    REQUIRE(
        TestSampleDistance(rng, inner_shell, neutron_group2) ==
        Approx(1. / 2.).epsilon(epsilon::Exponential(samples)));
  }
}
