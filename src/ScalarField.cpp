#include "ScalarField.hpp"

#include "Point.hpp"

#include <cassert>
#include <string>

// ScalarField

std::unique_ptr<const ScalarField>
ScalarField::Create(const pugi::xml_node& scalar_field_node) noexcept {
  std::unique_ptr<const ScalarField> scalar_field{}; // hope this triggers NRVO
  const std::string scalar_field_type = scalar_field_node.name();
  if (scalar_field_type == "constant"){
    scalar_field = std::make_unique<const ConstantField>(scalar_field_node);
  }
  else if (scalar_field_type == "linear"){
    scalar_field = std::make_unique<const LinearField>(scalar_field_node);
  }
  else {
    assert(false); // this should have been caught by the validator
  }
  return scalar_field;
}

// ConstantField

ConstantField::ConstantField(const pugi::xml_node& scalar_field_node) noexcept
    : c{scalar_field_node.attribute("constant").as_double()} {}

Real ConstantField::at(const Point&) const noexcept { return c; }

// LinearField

LinearField::LinearField(const pugi::xml_node& scalar_field_node) noexcept
    : g{scalar_field_node.child("gradient")},
      b{scalar_field_node.child("intercept").attribute("b").as_double()} {}

Real LinearField::at(const Point& p) const noexcept { return g.Dot(p) + b; }
