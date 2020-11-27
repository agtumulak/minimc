# https://github.com/catchorg/Catch2/blob/devel/docs/cmake-integration.md#cmake-target
Include(FetchContent)
FetchContent_Declare(
  catch_download
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.3
  )
FetchContent_MakeAvailable(catch_download)

# Add new unit test files here
set(UNIT_TEST_SOURCES
  test_hello_world.cpp
  )
# https://github.com/catchorg/Catch2/blob/devel/docs/slow-compiles.md
add_executable(all_tests test_main.cpp ${UNIT_TEST_SOURCES})
target_link_libraries(all_tests Catch2::Catch2 minimc)
list(APPEND CMAKE_MODULE_PATH ${catch_download_SOURCE_DIR}/contrib)
include(ParseAndAddCatchTests)
ParseAndAddCatchTests(all_tests)
