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
  /// @brief Virtual destructor (C++ Core Guidelines C.127)
  virtual ~ScalarField() noexcept;
  /// @brief Returns the value at a given Point
  virtual Real at(const Point& p) const noexcept = 0;
};

/// @brief A field where the value is independent of the Point
class ConstantField : public ScalarField {
public:
  /// @brief Constructs a constant field from a `constant` scalar field node
  ConstantField(const pugi::xml_node& scalar_field_node) noexcept;
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
  /// @brief Returns the linearly dependent value
  Real at(const Point& p) const noexcept override;
private:
  // Gradient of the field
  const Point g;
  // Value of the field at the origin
  const Real b;
};
