#include "Estimator.hpp"
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
  size_t samples = 1000; // set in point_source_leakage.xml
  const auto estimated_leakage_rate =
      f.Solve().FindEstimatorByName("leakage").GetScore(0, samples);
  auto expected_leakage_rate = std::exp(-1);
  REQUIRE_THAT(
      estimated_leakage_rate,
      Catch::Matchers::WithinRel(
          expected_leakage_rate,
          epsilon::Bernoulli(expected_leakage_rate, samples)));
}
