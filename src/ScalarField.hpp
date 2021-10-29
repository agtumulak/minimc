#pragma once

#include "BasicTypes.hpp"
#include "Point.hpp"

#include <memory>

/// @brief Abstract interface for all scalar fields
class ScalarField {
public:
  /// @brief Factory method to create new ScalarField from an XML document
  static std::unique_ptr<const ScalarField>
  Create(const pugi::xml_node& scalar_field_node) noexcept;
  /// @brief Constructs a ScalarField from by assigning member directly
  ScalarField(Real upper_bound, Real lower_bound) noexcept;
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~ScalarField() noexcept;
  /// @brief Returns true if the field is constant
  /// @details Used during input parsing to check if a TransportMethod is valid
  virtual bool IsConstant() const noexcept = 0;
  /// @brief Returns the value at a given Point
  virtual Real at(const Point& p) const noexcept = 0;
  /// @brief Upper bound on values that will be encountered (C++ Core
  ///        Guidelines C.131)
  const Real upper_bound;
  /// @brief Lower bound on values that will be encountered (C++ Core
  ///        Guidelines C.131)
  const Real lower_bound;
};

/// @brief A field where the value is independent of the Point
class ConstantField : public ScalarField {
public:
  /// @brief Constructs a constant field from a `constant` scalar field node
  ConstantField(const pugi::xml_node& scalar_field_node) noexcept;
  /// @brief Constructs a constant field from a given constant
  ConstantField(const Real c) noexcept;
  /// @brief Returns true because a ConstantField is constant
  bool IsConstant() const noexcept override;
  /// @brief Returns the constant value (C++ Core Guidelines F.9)
  Real at(const Point&) const noexcept override;

private:
  const Real c;
};

/// @brief A field where the value is only a linear function of Point
/// @details The value of the function @f$ f(\boldsymbol{p}) @f$ at some
///          point is just
///          @f[
///            f(\boldsymbol{p}) = \boldsymbol{g}^{T} \boldsymbol{p} + b
///          @f]
class LinearField : public ScalarField {
public:
  /// @brief Constructs a linear field from a `linear` scalar field node
  LinearField(const pugi::xml_node& scalar_field_node) noexcept;
  /// @brief Returns true iff the gradient of the linear field is zero
  bool IsConstant() const noexcept override;
  /// @brief Returns the linearly dependent value
  Real at(const Point& p) const noexcept override;
private:
  // Gradient of the field
  const Point g;
  // Value of the field at the origin
  const Real b;
};
