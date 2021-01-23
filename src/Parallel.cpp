#include "Parallel.hpp"

// ChunkGiver

//// public

ChunkGiver::ChunkGiver(size_t last, size_t chunksize)
    : last{last}, chunksize{chunksize} {}

std::optional<std::pair<size_t, size_t>> ChunkGiver::Next() {
  const std::lock_guard<std::mutex> lock(m);
  next_begin = next_end;
  if (next_begin > last) {
    return std::nullopt;
  }
  next_end = next_begin + chunksize;
  if (next_end > last + 1) {
    return std::make_pair(next_begin, last + 1);
  }
  return std::make_pair(next_begin, next_end);
}
