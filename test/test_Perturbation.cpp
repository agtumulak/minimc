#include "Estimator/Estimator.hpp"
#include "FixedSource.hpp"
#include "Perturbation/Sensitivity/Sensitivity.hpp"
#include "Statistics.hpp"
#include "World.hpp"
#include "XMLDocument.hpp"
#include "catch2/catch_test_macros.hpp"
#include "catch2/matchers/catch_matchers.hpp"
#include "catch2/matchers/catch_matchers_floating_point.hpp"

#include <cmath>

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

TEST_CASE(
    "nonmultiplying sphere leakage rate sensitivity to total cross section") {
  // compute sensitivity using differential operator sampling
  XMLDocument doc{"nonmultiplying_sphere_leakage_totalxs_sensitivity.xml"};
  const World w{doc.root};
  FixedSource f{doc.root};
  const auto dos =
      f.Solve().front()->sensitivities.front()->GetScore(0, f.total_weight);
  // perturb capture cross section by +0.05b
  XMLDocument doc_hi{"nonmultiplying_sphere_leakage_totalxs_perturb_hi.xml"};
  const World w_hi{doc_hi.root};
  FixedSource f_hi{doc_hi.root};
  const auto perturb_hi = f_hi.Solve().front()->GetScore(0, f_hi.total_weight);
  // perturb capture cross section by -0.05b
  XMLDocument doc_lo{"nonmultiplying_sphere_leakage_totalxs_perturb_lo.xml"};
  const World w_lo{doc_lo.root};
  FixedSource f_lo{doc_lo.root};
  const auto perturb_lo = f_lo.Solve().front()->GetScore(0, f_lo.total_weight);
  // compare differential operator sampling against central difference
  const auto central_difference = (perturb_hi - perturb_lo) / 0.1;
  REQUIRE_THAT(dos, Catch::Matchers::WithinRel(central_difference, 0.1));
}
