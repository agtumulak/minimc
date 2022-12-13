#include "Estimator/Estimator.hpp"
#include "FixedSource.hpp"
#include "Statistics.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cmath>
#include <cstddef>

TEST_CASE("point source leakage probability") {
  XMLDocument doc{"point_source_leakage.xml"};
  const World w{doc.root};
  FixedSource f{doc.root};
  const auto estimated_leakage_rate =
      f.Solve().front()->GetScore(0, f.total_weight);
  auto expected_leakage_rate = std::exp(-1);
  REQUIRE_THAT(
      estimated_leakage_rate,
      Catch::Matchers::WithinRel(
          expected_leakage_rate,
          epsilon::Bernoulli(expected_leakage_rate, f.total_weight)));
}
