#include "BasicTypes.hpp"
#include "Particle.hpp"
#include "Point.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_exception.hpp"

#include <variant>
#include <stdexcept>

TEST_CASE("Energy and Type constructor") {
  Particle p{Point{}, Direction{1, 0, 0}, Group{42}, Particle::Type::photon, 1};
  REQUIRE(p.GetPosition() == Point{0, 0, 0});
  REQUIRE(std::get<Group>(p.GetEnergy()) == Group{42});
  REQUIRE(p.type == Particle::Type::photon);
}

TEST_CASE("stream Particle some distance") {
  Particle p{Point{}, Direction{1, 0, 0}, Group{1}, Particle::Type::neutron, 1};
  p.Stream(10);
  REQUIRE(p.GetPosition() == Point{10, 0, 0});
}

TEST_CASE("Particle Cell can be assigned") {
  XMLDocument doc{"multigroup.xml"};
  const auto w{World{doc.root}};
  Particle p{Point{}, Direction{1, 0, 0}, Group{1}, Particle::Type::neutron, 1};
  REQUIRE_THROWS_MATCHES(
      p.GetCell(), std::runtime_error,
      Catch::Matchers::Message("Particle does not belong to a Cell"));
  p.SetCell(w.cells.front());
  REQUIRE(p.GetCell() == w.cells.front());
}
