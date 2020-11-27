#include "catch2/catch.hpp"
#include "hello_world.hpp"

TEST_CASE("simplest possible test case") { REQUIRE(true); }

TEST_CASE("HelloWorld class") {
  HelloWorld hello_world{};
  REQUIRE(hello_world.SayHello() == "Hello World!");
  REQUIRE(hello_world.SayGoodbye() == "Goodbye!");
}
