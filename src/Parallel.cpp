#include "Parallel.hpp"

// ChunkGiver

//// public

ChunkGiver::ChunkGiver(size_t size, size_t chunk) : size{size}, chunk{chunk} {}

std::optional<std::pair<size_t, size_t>> ChunkGiver::Next() {
  const std::lock_guard<std::mutex> lock(m);
  begin = end;
  if (begin >= size) {
    return std::nullopt;
  }
  end = begin + chunk;
  if (end > size) {
    return std::make_pair(begin, size);
  }
  return std::make_pair(begin, end);
}

void ChunkGiver::Reset(size_t new_size) noexcept {
  size = new_size;
  begin = 0;
  end = 0;
}
