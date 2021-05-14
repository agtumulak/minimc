#include "BasicTypes.hpp"
#include "Nuclide.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "Reaction.hpp"
#include "Statistics.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

#include <algorithm>
#include <cstddef>
#include <variant>
#include <vector>

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
  }
}
