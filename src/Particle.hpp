#pragma once

#include "Bank.hpp"
#include "BasicTypes.hpp"
#include "Point.hpp"
#include "Reaction.hpp"

#include <iosfwd>
#include <memory>
#include <vector>

namespace Perturbation {
namespace IndirectEffect {
class Interface;
}
} // namespace Perturbation

class Cell;

/// @brief The primary entity performing random walks in a World.
/// @details The awkward declaration order of member variables is meant to
///          improve alignment.
/// @note Users of this class should assume that Direction is normalized to
///       avoid extra computation.
/// @todo Make Particle a struct or make most members public. Change Particle
///       documentation from emphasis on encapsulation to emphasis on
///       const-correctness when references to it are used.
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
  /// @brief Member constructor. Explicitly assigns phase-space members.
  Particle(
      const std::vector<std::shared_ptr<
          Perturbation::IndirectEffect::Interface>>& indirect_effects,
      const Point& position, const Direction& direction, const Energy& energy,
      const Cell* cell, const Type type, RNG::result_type seed) noexcept;
  /// @brief Scatters the Particle with an outgoing direction and energy.
  /// @details Scattering is assumed to be azimuthally symmetric.
  /// @param mu The scattering cosine @f$ \mu @f$
  /// @param e The outgoing energy
  void Scatter(const Real& mu, const Energy& e) noexcept;
  /// @brief Return the current position of the Particle
  const Point& GetPosition() const noexcept;
  /// @brief Updates the current position of the Particle
  void SetPosition(const Point& p) noexcept;
  /// @brief Return the current direction of the Particle
  const Direction& GetDirection() const noexcept;
  /// @brief Set the current direction of the Particle
  void SetDirection(const Direction& d) noexcept;
  /// @brief Sets the direction to a random isotropic direction
  /// @note This should be replaced by a method which accepts scattering cosine
  void SetDirectionIsotropic() noexcept;
  /// @brief Returns the current energy of the Particle
  const Energy& GetEnergy() const noexcept;
  /// @brief Updates the current energy of the Particle
  void SetEnergy(const Energy& e) noexcept;
  /// @brief Returns the Type of the Particle
  Type GetType() const noexcept;
  /// @brief Returns a reference to the current Cell the Particle is within
  const Cell& GetCell() const;
  /// @brief Sets the current Cell occupied by the Particle
  void SetCell(const Cell& c) noexcept;
  /// @brief Banks secondaries produced during transport using an outgoing
  ///        Direction and outgoing Energy
  void
  BankSecondaries(const Direction& direction, const Energy& energy) noexcept;
  /// @brief Transfers secondaries produced by this Particle to the front of
  ///        a given Bank
  /// @details This is implemented using std::list::splice so this Particle's
  ///          own list of secondaries becomes empty after the operation.
  void MoveSecondariesTo(Bank& bank) noexcept;
  /// @brief Sample a random number uniformly in @f$ [0, 1) @f$
  /// @details This is provided as a convenience function because more
  ///          complicated sampling schemes require multiple uniformly
  ///          distributed random numbers
  Real Sample() noexcept;
  /// @brief Returns true if the Particle should continue to be transported
  bool IsAlive() const noexcept;
  /// @brief Indirect effects due to a Perturbation (C++ Core Guidelines C.131)
  const std::vector<std::shared_ptr<Perturbation::IndirectEffect::Interface>>
      indirect_effects;

private:
  // Secondaries produced
  Bank secondaries;
  // Position may be anywhere in @f$ \mathbb{R}^3 @f$
  Point position{0, 0, 0};
  // Direction must be constrained to @f$ \lVert v \rVert = 1 @f$
  Direction direction{1, 0, 0};
  // Energy in continuous energy calculation or group in multigroup calculation
  Energy energy{Group{1}};

public:
  /// @brief Non-owning pointer to current const Cell occupied by the Particle
  ///        (C++ Core Guideliens C.131)
  const Cell* cell{nullptr};
  /// @brief Type of the Particle (C++ Core Guidelines C.131)
  const Type type{Type::neutron};
  /// @brief Random number generator (C++ Core Guidelines C.131)
  /// @details This Particle contains its own random number generator. Any
  ///          Particle initialized with the same member variables and
  ///          transported in the same World should undergo the same history
  ///          and return the same Bank.
  RNG rng;
  /// @brief Flag describing the reaction that occured at the current point in
  ///        phase space
  Reaction reaction{Reaction::birth};
};
