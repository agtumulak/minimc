#include "Parallel.hpp"
#include "catch2/catch.hpp"

#include <cstddef>
#include <future>
#include <numeric>
#include <utility>
#include <vector>

TEST_CASE("split chunks") {
  // One chunk only
  ChunkGiver single_chunk{10, 1000};
  REQUIRE(single_chunk.Next() == std::pair{size_t{0}, size_t{10}});
  REQUIRE(single_chunk.Next() == std::nullopt);
  // perfectly even
  ChunkGiver perfectly_even{999, 333};
  REQUIRE(perfectly_even.Next() == std::pair{size_t{0}, size_t{333}});
  REQUIRE(perfectly_even.Next() == std::pair{size_t{333}, size_t{666}});
  REQUIRE(perfectly_even.Next() == std::pair{size_t{666}, size_t{999}});
  REQUIRE(perfectly_even.Next() == std::nullopt);
  // modulo one
  ChunkGiver mod_one{1000, 333};
  REQUIRE(mod_one.Next() == std::pair{size_t{0}, size_t{333}});
  REQUIRE(mod_one.Next() == std::pair{size_t{333}, size_t{666}});
  REQUIRE(mod_one.Next() == std::pair{size_t{666}, size_t{999}});
  REQUIRE(mod_one.Next() == std::pair{size_t{999}, size_t{1000}});
  REQUIRE(mod_one.Next() == std::nullopt);
  // modulo two
  ChunkGiver mod_two{1001, 333};
  REQUIRE(mod_two.Next() == std::pair{size_t{0}, size_t{333}});
  REQUIRE(mod_two.Next() == std::pair{size_t{333}, size_t{666}});
  REQUIRE(mod_two.Next() == std::pair{size_t{666}, size_t{999}});
  REQUIRE(mod_two.Next() == std::pair{size_t{999}, size_t{1001}});
  REQUIRE(mod_two.Next() == std::nullopt);
}

TEST_CASE("multithreaded chunk splitting") {
  size_t h{1000};
  size_t threads{10};
  // shared between threads
  ChunkGiver shared_chunk{h, 10};
  // expected result from each thread
  std::vector<std::future<size_t>> results;
  // spawn each thread asynchronously
  for (size_t i = 1; i <= threads; i++) {
    results.push_back(std::async(
        [](ChunkGiver& cg) {
          size_t accumulated = 0;
          while (auto range = cg.Next()) {
            const auto [begin, end] = range.value();
            accumulated += end - begin; // quick way of counting elements
          }
          return accumulated;
        },
        std::ref(shared_chunk)));
  }
  // accumulate asynchronously
  size_t sum = std::reduce(
      results.begin(), results.end(), size_t{0},
      [](const auto& accumulated, auto& future) {
        return accumulated + future.get();
      });
  REQUIRE(sum == h);
}
