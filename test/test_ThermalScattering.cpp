#include "FixedSource.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch_test_macros.hpp"

TEST_CASE("thermal scattering") {
  XMLDocument doc{"test_ThermalScattering.xml"};
  const World w{doc.root};
  FixedSource f{doc.root};
  // TODO: Check result of estimators, similar to `runsab` from
  //       `feature/sab-validation`
  REQUIRE_NOTHROW(f.Solve());
}
