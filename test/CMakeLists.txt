# https://github.com/catchorg/Catch2/blob/devel/docs/cmake-integration.md#cmake-target
Include(FetchContent)
FetchContent_Declare(
  catch_download
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.3
  )
FetchContent_MakeAvailable(catch_download)

# add new unit test source files here
set(UNIT_TEST_SOURCES
  test_CSGSurface.cpp
  test_Cell.cpp
  test_FixedSource.cpp
  test_HDF5.cpp
  test_Material.cpp
  test_Nuclide.cpp
  test_Particle.cpp
  test_Point.cpp
  test_ThermalScattering.cpp
  test_XMLDocument.cpp
  test_World.cpp
  )
# https://github.com/catchorg/Catch2/blob/devel/docs/slow-compiles.md
add_executable(all_tests catch_main.cpp ${UNIT_TEST_SOURCES})
target_link_libraries(all_tests Catch2::Catch2 minimc)

# add new unit test input files here
set(UNIT_TEST_INPUTS
  minimc.xsd
  minimc_basictypes.xsd
  point_source_leakage.xml
  simple_continuous.xml
  simple_multigroup.xml
  test_HDF5_valid.hdf5
  test_HDF5_invalid.hdf5
  test_XMLDocument.xsd
  test_XMLDocument_valid.xml
  test_XMLDocument_invalid.xml
  )
list(TRANSFORM UNIT_TEST_INPUTS PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/)
add_custom_target(
  copy_xml_inputs ALL
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${UNIT_TEST_INPUTS} ${CMAKE_CURRENT_BINARY_DIR}
  )

# replace CTest with Catch2
add_custom_target(test
  COMMAND all_tests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Run Catch2 instead of CTest")
add_dependencies(test all_tests copy_xml_inputs)

# code coverage with llvm
if (CMAKE_BUILD_TYPE MATCHES "Debug")
  add_custom_target(coverage
    COMMAND make test
    COMMAND xcrun llvm-profdata merge -sparse default.profraw -o default.profdata
    COMMAND xcrun llvm-cov show -instr-profile default.profdata all_tests -format=html -o report
    COMMAND open report/index.html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Create coverage report on macOS")
  add_dependencies(coverage all_tests)
endif()
