#include "BasicTypes.hpp"
#include "CSGSurface.hpp"
#include "Cell.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

TEST_CASE("Cell is constructed properly") {
  XMLDocument doc{"multigroup.xml"};
  const World w{doc.root};

  const auto& inner_shell = w.cells.at(1);
  const Particle neutron_group1{
      Point{}, Direction{1, 0, 0}, Group{1}, Particle::Type::neutron, 1};
  const Particle neutron_group2{
      Point{}, Direction{1, 0, 0}, Group{2}, Particle::Type::neutron, 1};

  SECTION("Cell returns correct nearest CSGSurface"){
    const auto& [surface, distance] =
        inner_shell.NearestSurface(Point{1.5, 0, 0}, Direction{1, 0, 0});
    REQUIRE(surface->name == "middle shell");
    REQUIRE_THAT(distance, Catch::Matchers::WithinRel(0.5));
  }
}
