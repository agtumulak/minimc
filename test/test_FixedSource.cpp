#include "BasicTypes.hpp"
#include "FixedSource.hpp"
#include "catch2/catch.hpp"

#include <future>
#include <numeric>
#include <vector>

using history = RNG::result_type;

TEST_CASE("split chunks") {
  // One chunk only
  ChunkGiver single_chunk{10, 1000};
  REQUIRE(single_chunk.Next() == std::pair{history{1}, history{11}});
  REQUIRE(single_chunk.Next() == std::nullopt);
  // perfectly even
  ChunkGiver perfectly_even{999, 333};
  REQUIRE(perfectly_even.Next() == std::pair{history{1}, history{334}});
  REQUIRE(perfectly_even.Next() == std::pair{history{334}, history{667}});
  REQUIRE(perfectly_even.Next() == std::pair{history{667}, history{1000}});
  REQUIRE(perfectly_even.Next() == std::nullopt);
  // modulo one
  ChunkGiver mod_one{1000, 333};
  REQUIRE(mod_one.Next() == std::pair{history{1}, history{334}});
  REQUIRE(mod_one.Next() == std::pair{history{334}, history{667}});
  REQUIRE(mod_one.Next() == std::pair{history{667}, history{1000}});
  REQUIRE(mod_one.Next() == std::pair{history{1000}, history{1001}});
  REQUIRE(mod_one.Next() == std::nullopt);
  // modulo two
  ChunkGiver mod_two{1001, 333};
  REQUIRE(mod_two.Next() == std::pair{history{1}, history{334}});
  REQUIRE(mod_two.Next() == std::pair{history{334}, history{667}});
  REQUIRE(mod_two.Next() == std::pair{history{667}, history{1000}});
  REQUIRE(mod_two.Next() == std::pair{history{1000}, history{1002}});
  REQUIRE(mod_two.Next() == std::nullopt);
}

TEST_CASE("multithreaded chunk splitting") {
  history h{1000};
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
