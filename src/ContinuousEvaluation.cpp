#include "ContinuousEvaluation.hpp"

#include "Constants.hpp"
#include "HDF5DataSet.hpp"
#include "pugixml.hpp"

// ContinuousEvaluation

//// public

ContinuousEvaluation::ContinuousEvaluation(
    const pugi::xml_node& evaluation_node)
    : xs{HDF5DataSet<1>{evaluation_node.attribute("file").as_string()}
             .ToContinuousMap()},
      temperature{evaluation_node.attribute("temperature").as_double()} {}

bool ContinuousEvaluation::IsValid(Temperature requested) const noexcept {
  const auto relative_error = (requested - temperature) / temperature;
  return relative_error < constants::relative_temperature_difference_tolerance;
}
