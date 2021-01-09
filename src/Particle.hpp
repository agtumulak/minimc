#pragma once

#include "BasicTypes.hpp"
#include "Point.hpp"

class Cell;

/// @brief The primary entity performing random walks in a World.
/// @details Particles are characterized by their position, direction, energy,
///          type, and an alive flag. The awkward declaration order of member
///          variables is meant to improve alignment. On Apple clang version
///          12.0.0 the members are as follows:
///          Member        | Size (bytes)  | Total (bytes)
///          ------------- | ------------- | -------------
///          position      |            24 |            24
///          direction     |            24 |            48
///          energy        |             8 |            56
///          cell*         |             8 |            64
///          type          |             4 |            68
///          alive         |             1 |            69
///          padding       |             3 |            72
/// @note All users of this class should assume that the direction member is
///       normalized to one to avoid extra computation. Maintainers of this
///       class must take care care to keep direction normalized.
class Particle {
public:
  /// @brief Affects which cross section data is used during transport, among
  ///        other things
  enum class Type {
    neutron,
    photon,
  };
  /// @brief Helper function to convert from std::string to Type
  static Type ToType(const std::string& name) noexcept;
  /// @brief Default constructor. Returns a neutron at the origin moving in
  ///        the x direction in Group 1.
  Particle() noexcept;
  /// @brief Type and Energy constructor. For accessing nuclear data.
  /// @details Similar to default constructor but with user-define
  ///          Particle::Type and Energy. These two parameters are all that
  ///          are used when accessing nuclear data.
  Particle(const Energy& energy, const Type type) noexcept;
  /// @brief Checks if `*this` is still alive
  explicit operator bool() const noexcept;
  /// @brief Return the current position of the Particle
  const Point& GetPosition() const noexcept;
  /// @brief Return the current direction of the Particle
  const Point& GetDirection() const noexcept;
  /// @brief Sets the direction to a random isotropic direction
  /// @note This should be replaced by a method which accepts scattering cosine
  void SetDirectionIsotropic(std::minstd_rand& rng) noexcept;
  /// @brief Returns the current energy of the Particle
  const Energy& GetEnergy() const noexcept;
  /// @brief Updates the current energy of the Particle
  void SetEnergy(const Energy& e) noexcept;
  /// @brief Returns a reference to the current Cell the Particle is within
  const Cell& GetCell() const;
  /// @brief Sets the current Cell occupied by the Particle
  void SetCell(const Cell& c) noexcept;

private:
  // Position may be anywhere in @f$ \mathbb{R}^3 @f$
  Point position{0, 0, 0};
  // Direction must be constrained to @f$ \lVert v \rVert = 1 @f$
  Point direction{1, 0, 0};
  // Energy in continuous energy calculation or group in multigroup calculation
  Energy energy{Group{1}};
  // Non-owning pointer to current const Cell occupied by the Particle
  const Cell* cell{nullptr};

public:
  /// @brief Type of the Particle (C++ Core Guidelines C.131)
  const Type type{Type::neutron};

private:
  // Flag for determining if this Particle is still alive
  bool alive{true};
};
