#include "BasicTypes.hpp"
#include "Nuclide.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_exception.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <stdexcept>

TEST_CASE("nonexistent nuclide name throws exception") {
  REQUIRE_THROWS_MATCHES(
      World{XMLDocument{"multigroup_nonexistent_nuclide.xml"}.root},
      std::runtime_error,
      Catch::Matchers::Message(
          "Nuclide node \"nonexistent\" not found. Must be one of: "
          "[\"nuclide A\", \"nuclide B\", ]"));
}

TEST_CASE("malformed multigroup data for Nuclide throws exception") {
  REQUIRE_THROWS_MATCHES(
      World{XMLDocument{"multigroup_malformed_capture.xml"}.root},
      std::runtime_error,
      Catch::Matchers::Message("/minimc/nuclides/multigroup/nuclide/neutron/"
                               "capture: Expected 2 entries "
                               "but got 3"));
  REQUIRE_THROWS_MATCHES(
      World{XMLDocument{"multigroup_malformed_scatter.xml"}.root},
      std::runtime_error,
      Catch::Matchers::Message("/minimc/nuclides/multigroup/nuclide/neutron/"
                               "scatter: Expected 4 entries "
                               "but got 5"));
  REQUIRE_THROWS_MATCHES(
      World{XMLDocument{"multigroup_missing_particle_data.xml"}.root},
      std::runtime_error,
      Catch::Matchers::Message(
          "/minimc/nuclides/multigroup/nuclide: \"neutron\" node not found"));
}

TEST_CASE("Nuclide member methods work properly") {
  SECTION("Multigroup methods") {
    XMLDocument doc{"multigroup.xml"};
    World w{doc.root};
    Particle neutron_group1{{},       Point{}, Direction{1, 0, 0},
                            Group{1}, {},      Particle::Type::neutron,
                            1};
    neutron_group1.SetCell(w.FindCellContaining(neutron_group1.GetPosition()));
    Particle neutron_group2{{},       Point{}, Direction{1, 0, 0},
                            Group{2}, {},      Particle::Type::neutron,
                            1};
    neutron_group2.SetCell(w.FindCellContaining(neutron_group2.GetPosition()));
    const Nuclide hydrogen{
        doc.root
            .select_node(
                "/minimc/nuclides/multigroup/nuclide[@name='hydrogen']")
            .node()};
    const Nuclide oxygen{
        doc.root
            .select_node("/minimc/nuclides/multigroup/nuclide[@name='oxygen']")
            .node()};
    const Nuclide uranium235{
        doc.root
            .select_node(
                "/minimc/nuclides/multigroup/nuclide[@name='uranium235']")
            .node()};

    REQUIRE(hydrogen.GetCellMajorant(neutron_group1) == 1);
    REQUIRE(hydrogen.GetCellMajorant(neutron_group2) == 1);
    REQUIRE(oxygen.GetCellMajorant(neutron_group1) == 1.5);
    REQUIRE(oxygen.GetCellMajorant(neutron_group2) == 1);
    REQUIRE(uranium235.GetCellMajorant(neutron_group1) == 1.33);
    REQUIRE(uranium235.GetCellMajorant(neutron_group2) == 2.67);
  }
  SECTION("Continuous methods") {
    XMLDocument doc{"continuous.xml"};
    World w{doc.root};
    Particle neutron{
        {},
        Point{},
        Direction{1, 0, 0},
        ContinuousEnergy{0.999e-6},
        {},
        Particle::Type::neutron,
        1};
    neutron.SetCell(w.FindCellContaining(neutron.GetPosition()));
    const Nuclide hydrogen{
        doc.root
            .select_node(
                "/minimc/nuclides/continuous/nuclide[@name='hydrogen']")
            .node()};
    const Nuclide oxygen{
        doc.root
            .select_node("/minimc/nuclides/continuous/nuclide[@name='oxygen']")
            .node()};
    const Nuclide uranium235{
        doc.root
            .select_node(
                "/minimc/nuclides/continuous/nuclide[@name='uranium235']")
            .node()};
    REQUIRE_THROWS_MATCHES(
        World{XMLDocument{"continuous_invalid_nuclide_data_path.xml"}.root},
        std::runtime_error,
        Catch::Matchers::Message("File not found: /not/real"));

    REQUIRE_THAT(
        hydrogen.GetCellMajorant(neutron),
        Catch::Matchers::WithinRel(20.74794400216));
    REQUIRE_THAT(
        oxygen.GetCellMajorant(neutron),
        Catch::Matchers::WithinRel(3.7969596274296));
    REQUIRE_THAT(
        uranium235.GetCellMajorant(neutron),
        Catch::Matchers::WithinRel(92.3574585600001));
  }
}
