#include "Estimator.hpp"

#include <ostream>

// Estimator

std::ostream& operator<<(std::ostream& os, const Estimator& e) noexcept {
  for (const auto& [event, score] : e) {
    os << Estimator::ToString(event) + ": " + std::to_string(score)
       << std::endl;
  }
  return os;
}

//// public

std::string Estimator::ToString(const Event e) noexcept {
  switch (e) {
  case Estimator::Event::capture:
    return "capture";
  case Estimator::Event::collision:
    return "collision";
  case Estimator::Event::fission:
    return "fission";
  case Estimator::Event::implicit_fission:
    return "implicit fission";
  case Estimator::Event::scatter:
    return "scatter";
  case Estimator::Event::surface_crossing:
    return "surface crossing";
  }
}

Real& Estimator::at(Event e) { return elements.at(e); }

const Real& Estimator::at(Event e) const { return elements.at(e); }

Estimator& Estimator::operator+=(const Estimator& rhs) noexcept {
  for (const auto& [event, score] : rhs) {
    elements.at(event) += score;
  }
  return *this;
}

Estimator::elements_type::const_iterator Estimator::begin() const noexcept {
  return elements.cbegin();
}
Estimator::elements_type::const_iterator Estimator::end() const noexcept {
  return elements.cend();
}
