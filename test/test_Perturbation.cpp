#include "Estimator/Estimator.hpp"
#include "Perturbation/Sensitivity/Sensitivity.hpp"
#include "FixedSource.hpp"
#include "Statistics.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cmath>
#include <cstddef>

TEST_CASE(
    "purely absorbing sphere leakage rate sensitivity to total cross section") {
  XMLDocument doc{"absorbing_sphere_leakage_totalxs_sensitivity.xml"};
  const World w{doc.root};
  FixedSource f{doc.root};
  const auto estimated_leakage_sensitivity =
      f.Solve().front()->sensitivities.front()->GetScore(0, f.total_weight);
  auto expected_sensitivity = -std::exp(-1);
  REQUIRE_THAT(
      estimated_leakage_sensitivity,
      Catch::Matchers::WithinRel(
          expected_sensitivity,
          epsilon::Bernoulli(-expected_sensitivity, f.total_weight)));
}
