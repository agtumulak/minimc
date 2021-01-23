#pragma once

#include <cstddef>
#include <mutex>
#include <optional>
#include <utility>

/// @brief Splits a range of nonnegative integers into fixed-size chunks.
/// @details The last chunk may have fewer values than `chunk`. Each call
///          to Next() will return a distinct subrange in a thread-safe manner.
class ChunkGiver {
public:
  /// @brief Constructs a ChunkGiver for `size` elements and fixed chunk size
  ChunkGiver(size_t size, size_t chunk);
  /// @brief Disallow copy/move constructor:
  ///        https://www.stroustrup.com/C++11FAQ.html#default2
  ChunkGiver(const ChunkGiver&) = delete;
  /// @brief Disallow copy/move assignment operator:
  ///        https://www.stroustrup.com/C++11FAQ.html#default2
  ChunkGiver& operator=(const ChunkGiver&) = delete;
  /// @brief Returns the next chunk. Returns std::nullopt if there are zero
  ///        chunks left to give.
  /// @return A pair of `[start, end)` integers representing a subrange.
  std::optional<std::pair<size_t, size_t>> Next();

private:
  std::mutex m;
  const size_t size;
  const size_t chunk;
  size_t begin{0};
  size_t end{0};
};
