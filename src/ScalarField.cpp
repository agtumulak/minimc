#include "ScalarField.hpp"

#include "Point.hpp"
#include "pugixml.hpp"

#include <cassert>
#include <iosfwd>
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

ScalarField::ScalarField(Real upper_bound, Real lower_bound) noexcept
    : upper_bound{upper_bound}, lower_bound{lower_bound} {}

ScalarField::~ScalarField() noexcept {}

// ConstantField

ConstantField::ConstantField(const pugi::xml_node& scalar_field_node) noexcept
    : ScalarField{
      scalar_field_node.attribute("c").as_double(),
      scalar_field_node.attribute("c").as_double()},
      c{scalar_field_node.attribute("c").as_double()} {}

ConstantField::ConstantField(const Real c) noexcept : ScalarField{c, c}, c{c} {}

bool ConstantField::IsConstant() const noexcept { return true; }

Real ConstantField::at(const Point&) const noexcept { return c; }

// LinearField

LinearField::LinearField(const pugi::xml_node& scalar_field_node) noexcept
    : ScalarField{
      scalar_field_node.child("bounds").attribute("upper").as_double(),
      scalar_field_node.child("bounds").attribute("lower").as_double()},
      g{scalar_field_node.child("gradient")},
      b{scalar_field_node.child("intercept").attribute("b").as_double()} {}

bool LinearField::IsConstant() const noexcept { return g == Point{0, 0, 0}; }

Real LinearField::at(const Point& p) const noexcept { return g.Dot(p) + b; }
