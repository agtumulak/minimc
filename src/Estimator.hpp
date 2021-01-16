#pragma once

#include <iostream>
#include <map>
#include <string>

/// @brief Keeps track of milestones in a Particle object's random walk
class Estimator {
  /// @brief Writes contents of Estimator to ostream
  friend std::ostream&
  operator<<(std::ostream& os, const Estimator& e) noexcept;

public:
  /// @brief Label for significant events
  enum class Event {
    surface_crossing,
    collision,
  };
  /// @brief Helper function to convert from Event to std::string
  static std::string ToString(const Event e) noexcept;

private:
  using elements_type = std::map<Event, size_t>;

public:
  /// @brief Returns a reference to a given Event
  size_t& at(Event e);
  /// @brief Returns a const reference to a given Event
  const size_t& at(Event e) const;
  /// @brief Adds the contents of rhs to the Estimator
  /// @details Only keys which exist in this Estimator are added to
  Estimator& operator+=(const Estimator& rhs);

private:
  // Used by ranged-based for loops
  elements_type::const_iterator begin() const noexcept;
  // Used by ranged-based for loops
  elements_type::const_iterator end() const noexcept;
  // Default initialized to zero counts
  elements_type elements{{Event::surface_crossing, 0}, {Event::collision, 0}};
};
