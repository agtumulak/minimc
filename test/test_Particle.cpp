#include "Particle.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

TEST_CASE("Particle default constructor") {
  Particle p{};
  REQUIRE(p); // checks operator bool()
  REQUIRE(p.GetPosition() == Point{0, 0, 0});
  REQUIRE(p.GetDirection() == Point{1, 0, 0});
}

TEST_CASE("Particle Cell can be assigned"){
  XMLDocument doc{"simple.xml"};
  const auto w{World{doc.root}};
  Particle p{};
  REQUIRE_THROWS_WITH(p.GetCell(), "Particle does not belong to a Cell");
  p.SetCell(w.cells.front());
  REQUIRE(p.GetCell() == w.cells.front());
}
