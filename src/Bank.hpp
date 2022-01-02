#pragma once

#include <list>
#include <type_traits>

class Particle;

/// @brief A collection of Particle objects
class Bank {
public:
  /// @brief Constant-time concatenation
  Bank& operator+=(Bank&& rhs) noexcept;
  /// @brief Returns a reference to the last Particle
  Particle& back() noexcept;
  /// @brief Returns true if the Bank has no Particle elements
  bool empty() const noexcept;
  /// @brief Constructs a Particle in-place
  template <typename... Args> Particle& emplace_back(Args&&... args) noexcept {
    return banked.emplace_back(std::forward<Args>(args)...);
  }
  /// @brief Removes the last Particle
  void pop_back() noexcept;

private:
  std::list<Particle> banked;
};
