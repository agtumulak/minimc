# https://github.com/catchorg/Catch2/blob/devel/docs/cmake-integration.md#cmake-target
Include(FetchContent)
FetchContent_Declare(
  catch_download
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.3
  )
FetchContent_MakeAvailable(catch_download)

# Add new unit test source files here
set(UNIT_TEST_SOURCES
  test_ValidateXML.cpp
  )
# https://github.com/catchorg/Catch2/blob/devel/docs/slow-compiles.md
add_executable(all_tests catch_main.cpp ${UNIT_TEST_SOURCES})
target_link_libraries(all_tests Catch2::Catch2 minimc)

# Add new unit test input files here
set(UNIT_TEST_INPUTS
  minimc_schema.xsd
  test_ValidateXML_valid.xml
  test_ValidateXML_invalid.xml
  )
list(TRANSFORM UNIT_TEST_INPUTS PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)
add_custom_target(
  copy_xml_inputs ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${UNIT_TEST_INPUTS} ${CMAKE_CURRENT_BINARY_DIR}
  )

# Replace CTest with Catch2
add_custom_target(test
  COMMAND all_tests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Run Catch2 instead of CTest")
add_dependencies(test all_tests)
