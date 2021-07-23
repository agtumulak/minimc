#include "ThermalScattering.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

TEST_CASE("create TSL data") {
  XMLDocument doc{"simple_continuous.xml"};
  const auto tsl_node =
      doc.root
          .select_node("nuclides/continuous/nuclide[@name='hydrogen']/neutron/"
                       "scatter/tsl")
          .node();
  REQUIRE_NOTHROW(ThermalScattering{tsl_node});
}
