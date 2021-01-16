#pragma once

#include "BasicTypes.hpp"

#include <mutex>
#include <optional>
#include <utility>

namespace parallel {

/// @brief Splits histories into fixed-size chunks.
/// @details The last chunk may have fewer values than `chunksize`. Each call
///          to Next() will return a distinct range in a thread-safe manner.
class ChunkGiver {
public:
  /// @brief Constructs a ChunkGiver with a last History and fixed chunk size
  ChunkGiver(History last, size_t chunksize);
  /// @brief Disallow copy/move constructor:
  ///        https://www.stroustrup.com/C++11FAQ.html#default2
  ChunkGiver(const ChunkGiver&) = delete;
  /// @brief Disallow copy/move assignment operator:
  ///        https://www.stroustrup.com/C++11FAQ.html#default2
  ChunkGiver& operator=(const ChunkGiver&) = delete;
  /// @brief Returns the next chunk. Returns std::nullopt if there are zero
  ///        chunks left to give.
  /// @return A pair of `[start, end)` History numbers.
  std::optional<std::pair<History, History>> Next();

private:
  std::mutex m;
  const History last;
  const size_t chunksize;
  History next_begin{1};
  History next_end{1};
};

} // namespace parallel
