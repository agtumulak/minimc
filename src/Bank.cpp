#include "Bank.hpp"

Bank& Bank::operator+=(Bank&& rhs) noexcept {
  banked.splice(banked.begin(), rhs.banked);
  return *this;
}

Particle& Bank::back() noexcept { return banked.back(); }

bool Bank::empty() const noexcept { return banked.empty(); }
