#include "Particle.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

#include <variant>

TEST_CASE("Particle default constructor") {
  Particle p{};
  REQUIRE(p.IsAlive()); // checks operator bool()
  REQUIRE(p.GetPosition() == Point{0, 0, 0});
  REQUIRE(p.GetDirection() == Point{1, 0, 0});
  REQUIRE(std::get<Group>(p.GetEnergy()) == Group{1});
  REQUIRE(p.type == Particle::Type::neutron);
}

TEST_CASE("Energy and Type constructor") {
  Particle p{Group{42}, Particle::Type::photon};
  REQUIRE(std::get<Group>(p.GetEnergy()) == Group{42});
  REQUIRE(p.type == Particle::Type::photon);
}

TEST_CASE("stream Particle some distance") {
  Particle p{};
  p.Stream(10);
  REQUIRE(p.GetPosition() == Point{10, 0, 0});
}

TEST_CASE("Particle can be killed") {
  Particle p{};
  REQUIRE(p.IsAlive() == true);
  p.Kill();
  REQUIRE(p.IsAlive() == false);
}

TEST_CASE("Particle Cell can be assigned") {
  XMLDocument doc{"simple.xml"};
  const auto w{World{doc.root}};
  Particle p{};
  REQUIRE_THROWS_WITH(p.GetCell(), "Particle does not belong to a Cell");
  p.SetCell(w.cells.front());
  REQUIRE(p.GetCell() == w.cells.front());
}
