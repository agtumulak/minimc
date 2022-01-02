#pragma once

#include "Bank.hpp"
#include "BasicTypes.hpp"
#include "Point.hpp"

#include <iosfwd>
#include <map>
#include <memory>

class Cell;
class CSGSurface;
class EstimatorSetProxy;
class Nuclide;
class Perturbation;
class PerturbationSet;
class TransportMethod;
class World;

/// @brief The primary entity performing random walks in a World.
/// @details Particles are characterized by their position, direction, energy,
///          type, and an alive flag. The awkward declaration order of member
///          variables is meant to improve alignment.
/// @note Users of this class should assume that Direction is normalized to
///       avoid extra computation.
class Particle {
public:
  /// @brief Affects which cross section data is used during transport, among
  ///        other things
  enum class Type {
    neutron,
    photon,
  };
  /// @brief Mutually exclusive events which can occur
  enum class Event {
    birth,
    scatter,
    capture,
    fission,
    surface_cross,
    leak,
    virtual_collision,
  };
  /// @brief Helper function to convert from std::string to Type
  static Type ToType(const std::string& name) noexcept;
  /// @brief The TransportMethod used by all Particle objects
  static std::unique_ptr<const TransportMethod> transport_method;
  /// @brief Member constructor. Explicitly assigns phase-space members.
  Particle(
      const Point& position, const Direction& direction, const Energy& energy,
      const Type type, RNG::result_type seed,
      const Cell* cell = nullptr) noexcept;
  /// @brief Use Particle::transport_method to update Particle state until it
  ///        dies
  Bank Transport(EstimatorSetProxy& e, const World& w) noexcept;
  /// @brief Moves the particle along its current direction a given distance
  void Stream(const Real distance) noexcept;
  /// @brief Scatters the Particle with an outgoing direction and energy.
  /// @details Scattering is assumed to be azimuthally symmetric.
  /// @param mu The scattering cosine @f$ \mu @f$
  /// @param e The outgoing energy
  void Scatter(const Real& mu, const Energy& e) noexcept;
  /// @brief Returns a reference to the indirect effects of this Particle
  Real GetIndirectEffect(const Perturbation* perturbation) const noexcept;
  /// @brief Sets the indirect effects of the Particle
  void SetPerturbations(const PerturbationSet& perturbations) noexcept;
  /// @brief Return the current position of the Particle
  const Point& GetPosition() const noexcept;
  /// @brief Return the current direction of the Particle
  const Direction& GetDirection() const noexcept;
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
  /// @brief Sample a Nuclide given that the Particle has collided inside its
  ///        Cell
  const Nuclide& SampleNuclide() noexcept;
  /// @brief Returns true if the Particle should continue to be transported
  bool IsAlive() const noexcept;

private:
  // Indirect effects due to a Perturbation
  std::map<const Perturbation*, Real> indirect_effects;
  // Secondaries produced
  Bank secondaries;
  // Position may be anywhere in @f$ \mathbb{R}^3 @f$
  Point position{0, 0, 0};
  // Direction must be constrained to @f$ \lVert v \rVert = 1 @f$
  Direction direction{1, 0, 0};
  // Energy in continuous energy calculation or group in multigroup calculation
  Energy energy{Group{1}};
  // Non-owning pointer to current const Cell occupied by the Particle
  const Cell* cell{nullptr};

public:
  /// @brief Random number generator (C++ Core Guidelines C.131)
  /// @details This Particle contains its own random number generator. Any
  ///          Particle initialized with the same member variables and passed
  ///          to the same TransportMethod::Transport call should return the
  ///          same TransportMethod::Outcome.
  RNG rng;
  /// @brief Type of the Particle (C++ Core Guidelines C.131)
  const Type type{Type::neutron};
  /// @brief Flag describing the event that occured at the current point in
  ///        phase space
  Event event{Event::birth};
  /// @brief Pointer to most recent surface crossed
  std::shared_ptr<const CSGSurface> current_surface;
};
