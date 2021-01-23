#pragma once

#include <cstddef>
#include <mutex>
#include <optional>
#include <utility>

/// @brief Splits a contiguous range of integers into fixed-size chunks.
/// @details The last chunk may have fewer values than `chunksize`. Each call
///          to Next() will return a distinct subrange in a thread-safe manner.
class ChunkGiver {
public:
  /// @brief Constructs a ChunkGiver with a last integer and fixed chunk size
  ChunkGiver(size_t last, size_t chunksize);
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
  const size_t last;
  const size_t chunksize;
  size_t next_begin{1};
  size_t next_end{1};
};
