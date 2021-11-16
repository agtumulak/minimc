#include "History.hpp"

// History

//// public

History::History(
    const World& world, const Source& source, RNG::result_type seed) noexcept {
  history.emplace_back(world, source, seed);
}

void History::CrossSurface(
    std::shared_ptr<const CSGSurface> surface, const Real distance,
    const Cell& cell) noexcept {
  history.emplace_back(GetState(), surface, distance, cell);
}

void History::CollideWithinCell(const Real distance) noexcept {
  history.emplace_back(GetState(), distance);
}

void History::StreamWithinCell(const Real distance) noexcept {
  history.back().position += history.back().direction * distance;
}

const State& History::GetState() const noexcept { return history.back(); }

RNG& History::GetRNG() noexcept { return history.back().rng; }

std::vector<State>::const_iterator History::begin() const noexcept {
  return history.cbegin();
}

std::vector<State>::const_iterator History::end() const noexcept {
  return history.cend();
}
