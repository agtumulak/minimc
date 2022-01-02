#include "Estimator.hpp"
#include "FixedSource.hpp"
#include "Statistics.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch.hpp"

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
  REQUIRE(
      estimated_leakage_rate ==
      Approx(expected_leakage_rate)
          .epsilon(epsilon::Bernoulli(expected_leakage_rate, samples)));
}
