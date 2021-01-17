#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

#include <algorithm>

TEST_CASE("construct a World") {
  XMLDocument doc{"simple_multigroup.xml"};
  REQUIRE_NOTHROW(World{doc.root});
}

TEST_CASE("World is constructed properly") {
  // Test interdependencies between objects in a World (Cell on Material,
  // Material on Nuclide, Cell on CSGSurface). Tests on standalone objects
  // (CSGsurface, Nuclide) or stateless static members are performed elsewhere.
  XMLDocument doc{"simple_multigroup.xml"};
  const World w{doc.root};

  // only three Cells were declared in the input file
  REQUIRE(w.cells.size() == 3);
  const auto& pit = w.cells.at(0);
  const auto& inner_shell = w.cells.at(1);
  const auto& outer_shell = w.cells.at(2);
  REQUIRE(pit.name == "pit");
  REQUIRE(inner_shell.name == "inner shell");
  REQUIRE(outer_shell.name == "outer shell");

  // test dependence of multiple Cell objects on the same Material
  REQUIRE(pit.material == outer_shell.material);
  REQUIRE(pit.material != inner_shell.material);

  // test dependence of multiple Cell objects on the same CSGSurface
  const auto& pit_outer_sphere = *std::find_if(
      pit.surface_senses.cbegin(), pit.surface_senses.cend(),
      [](const auto& pair) { return pair.first->name == "inner sphere"; });
  const auto& inner_shell_inner_sphere = *std::find_if(
      inner_shell.surface_senses.cbegin(), inner_shell.surface_senses.cend(),
      [](const auto& pair) { return pair.first->name == "inner sphere"; });
  REQUIRE(pit_outer_sphere.first == inner_shell_inner_sphere.first);

  SECTION("FindCellContaining() finds correct Cell containing Point"){
    REQUIRE(w.FindCellContaining(Point{0,0,0}) == pit);
    REQUIRE(w.FindCellContaining(Point{0.5,0,0}) == pit);
    REQUIRE(w.FindCellContaining(Point{1,0,0}) == inner_shell);
    REQUIRE(w.FindCellContaining(Point{1.5,0,0}) == inner_shell);
    REQUIRE(w.FindCellContaining(Point{2,0,0}) == outer_shell);
    REQUIRE(w.FindCellContaining(Point{2.5,0,0}) == outer_shell);
  }
}
